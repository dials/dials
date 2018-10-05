#!/usr/bin/env python
# coding: utf-8
from __future__ import absolute_import, division, print_function

help_message =  """
This program performs scaling on integrated datasets, which attempts to improve
the internal consistency of the reflection intensities by correcting for
various experimental effects. By default, a physical scaling model is used,
with scale, decay and absorption components. If multiple input files have been
specified, the datasets will be jointly scaled against a common target of
unique reflection intensities.

The program outputs one scaled.pickle and scaled_experiments.json file, which
contains reflection data and scale models, from one or more experiments.
The output pickle file contains intensity.scale.value, the unscaled intensity
values used to determine the scaling model, and a inverse scale factor per
reflection. These values can then be used to merge the data for downstream
structural solution. Alternatively, the scaled_experiments.json and
scaled.pickle files can be passed back to dials.scale, and further scaling will
be performed, starting from where the previous job finished.

The scaling models determined by this program can be plotted with::

  dials.plot_scaling_models scaled.pickle scaled_experiments.json

Example use cases

Regular single-sweep scaling, with no absorption correction::

  dials.scale integrated.pickle integrated_experiments.json absorption_term=False

Scaling multiple datasets, specifying scale parameter interval::

  dials.scale 1_integrated.pickle 1_integrated_experiments.json 2_integrated.pickle 2_integrated_experiments.json scale_interval=10.0

Incremental scaling (with different options per dataset)::

  dials.scale integrated.pickle integrated_experiments.json scale_interval=10.0

  dials.scale integrated_2.pickle integrated_experiments_2.json scaled.pickle scaled_experiments.json scale_interval=15.0

"""
import time
import logging
import sys
import gc
from libtbx import phil
from libtbx.utils import Sorry
from libtbx.str_utils import make_sub_header
from cctbx import miller, crystal
import iotbx.merging_statistics
from dials.util import halraiser, log
from dials.array_family import flex
from dials.util.options import OptionParser, flatten_reflections,\
  flatten_experiments
from dials.util.version import dials_version
from dials.algorithms.scaling.scaling_library import create_scaling_model,\
  create_datastructures_for_structural_model, create_datastructures_for_target_mtz
from dials.algorithms.scaling.scaler_factory import create_scaler,\
  MultiScalerFactory
from dials.algorithms.scaling.scaling_utilities import parse_multiple_datasets,\
  select_datasets_on_ids, save_experiments, save_reflections,\
  assign_unique_identifiers, log_memory_usage, DialsMergingStatisticsError
from dials.algorithms.scaling.post_scaling_analysis import \
  exclude_on_batch_rmerge, exclude_on_image_scale


logger = logging.getLogger('dials')
info_handle = log.info_handle(logger)
phil_scope = phil.parse('''
  debug = False
    .type = bool
    .help = "Output additional debugging information"
  model = *physical array KB
      .type = choice
      .help = "Set scaling model to be applied to input datasets without
               an existing model. "
  output {
    log = dials.scale.log
      .type = str
      .help = "The log filename"
    debug.log = dials.scale.debug.log
      .type = str
      .help = "The debug log filename"
    plot_scaling_models = False
      .type = bool
      .help = "Option to switch on plotting of the scaling models determined."
    experiments = "scaled_experiments.json"
      .type = str
      .help = "Option to set filepath for output json."
    reflections = "scaled.pickle"
      .type = str
      .help = "Option to set filepath for output pickle file of scaled
               intensities."
    unmerged_mtz = None
      .type = path
      .help = "Filename to export an unmerged_mtz, calls dials.export internally."
    merged_mtz = None
      .type = path
      .help = "Filename to export a merged_mtz file."
    crystal_name = XTAL
      .type = str
      .help = "The crystal name to be exported in the mtz file metadata"
    use_internal_variance = False
      .type = bool
      .help = "Option to use internal spread of the intensities when merging
              reflection groups and calculating sigI, rather than using the
              sigmas of the individual reflections."
    merging.nbins = 20
      .type = int
      .help = "Number of bins to use for calculating and plotting merging stats."
    exclude_on_image_scale = None
      .type = float
      .help = "If set, images where the image inverse scale (defined by the
              scale component of the relevant model) is below this
              value will be set as outliers, and not included in merging stats
              or output for downstream processing. This is performed after the
              scaling algorithm has been run in the 'post-scaling' step.
              This option is intended as a quick way to exclude radiation
              damaged images, based on the principle that images that need
              to be significantly scaled up (relative to others) are likely to
              have unreliable intensities and scales. This option should be most
              appropriate for scaling multi-dataset thin-wedge datasets that
              are expected to have significant radiation damage."
    exclude_on_batch_rmerge = None
      .type = float
      .help = "If set, images which have an Rmerge above this value will be set
              as outliers, and not included in merging stats or output for
              downstream processing. This is performed after the scaling
              algorithm has been run in the 'post-scaling' step."
  }
  include scope dials.algorithms.scaling.scaling_options.phil_scope
  include scope dials.algorithms.scaling.cross_validation.cross_validate.phil_scope
  include scope dials.algorithms.scaling.scaling_refiner.scaling_refinery_phil_scope
''', process_includes=True)


class Script(object):
  """Main script to run the scaling algorithm."""

  def __init__(self, params, experiments, reflections):
    self.params = params
    self.experiments = experiments
    self.reflections = reflections
    self.scaler = None
    self.scaled_miller_array = None
    self.dataset_ids = []
    self.merging_statistics_result = None
    logger.debug('Initialised scaling script object')
    log_memory_usage()

  def run(self, save_data=True):
    """Run the scaling script."""
    start_time = time.time()
    self.prepare_input()
    self.scale()
    try:
      self.merging_stats()
    except DialsMergingStatisticsError as e:
      logger.info(e)
    if save_data:
      self.output()
    # All done!
    finish_time = time.time()
    logger.info("\nTotal time taken: {0:.4f}s ".format(finish_time - start_time))
    logger.info('\n'+'='*80+'\n')

  def prepare_input(self):
    """Perform checks on the data and prepare the data for scaling."""
    if len(self.experiments) != 1:
      logger.info('Checking for the existence of a reflection table \n'
        'containing multiple scaled datasets \n')
      self.reflections, self.dataset_ids = parse_multiple_datasets(self.reflections)
      logger.info("Found %s reflection tables in total.", len(self.reflections))
      logger.info("Found %s experiments in total.", len(self.experiments))

    if (self.params.dataset_selection.use_datasets or
     self.params.dataset_selection.exclude_datasets):
     if self.experiments.identifiers().count('') > 0:
      logger.warn('\nERROR: Attempting to choose datasets based on unique identifier,\n'
        'but not all datasets currently have a unique identifier! Please make\n'
        'sure all identifiers are set before attempting to select datasets.\n')
      logger.info('Current identifiers set as: %s', list(
        self.experiments.identifiers()))
      sys.exit()

    self.experiments, self.reflections = assign_unique_identifiers(
      self.experiments, self.reflections)
    logger.info("\nDataset unique identifiers are %s \n", list(
      self.experiments.identifiers()))

    if (self.params.dataset_selection.use_datasets or
     self.params.dataset_selection.exclude_datasets):
      self.experiments, self.reflections = \
        select_datasets_on_ids(self.experiments, self.reflections,
          use_datasets=self.params.dataset_selection.use_datasets,
          exclude_datasets=self.params.dataset_selection.exclude_datasets)

      logger.info("\nDataset unique identifiers for retained datasets are %s \n",
        list(self.experiments.identifiers()))


    if self.params.scaling_options.space_group:
      for experiment in self.experiments:
        sg_from_file = experiment.crystal.get_space_group().info()
        s_g_symbol = self.params.scaling_options.space_group
        crystal_symmetry = crystal.symmetry(space_group_symbol=s_g_symbol)
        experiment.crystal.set_space_group(crystal_symmetry.space_group())
        if crystal_symmetry.space_group() != sg_from_file:
          msg = ('WARNING: Manually overriding space group from {0} to {1}. {sep}'
            'If the reflection indexing in these space groups is different, {sep}'
            'bad things may happen!!! {sep}').format(sg_from_file, s_g_symbol, sep='\n')
          logger.info(msg)
    s_g_1 = self.experiments[0].crystal.get_space_group()
    for experiment in self.experiments:
      if experiment.crystal.get_space_group().type().number() \
        != s_g_1.type().number():
        logger.info("Space groups not determined to be equal; %s and %s", s_g_1.info(),
          experiment.crystal.get_space_group().info())
        raise Sorry('experiments have different space groups and cannot be '
          'scaled together, please reanalyse the data so that the space groups '
          'are consistent or manually specify a space group. Alternatively, '
          'some datasets can be excluded using the option exclude_datasets=')

    logger.info("Space group being used during scaling is %s" % s_g_1.info())

    if self.params.scaling_options.target_model:
      logger.info("Extracting data from structural model.")
      exp, reflections = create_datastructures_for_structural_model(
        self.reflections, self.experiments, self.params.scaling_options.target_model)
      self.experiments.append(exp)
      self.reflections.append(reflections)
      self.reflections[-1]['id'] = flex.int(self.reflections[-1].size(),
        len(self.reflections))

    elif self.params.scaling_options.target_mtz:
      logger.info("Extracting data from merged mtz.")
      exp, reflections = create_datastructures_for_target_mtz(self.experiments,
        self.params.scaling_options.target_mtz)
      self.experiments.append(exp)
      self.reflections.append(reflections)
      self.reflections[-1]['id'] = flex.int(self.reflections[-1].size(),
        len(self.reflections))

    if len(self.experiments) != len(self.reflections):
      raise Sorry("Mismatched number of experiments and reflection tables found.")

    # Perform any cutting of the dataset prior to creating scaling models.
    for reflection in self.reflections:
      reflection.set_flags(flex.bool(reflection.size(), False),
        reflection.flags.user_excluded_in_scaling)
      if self.params.cut_data.d_min:
        reflection.set_flags(reflection['d'] < self.params.cut_data.d_min,
        reflection.flags.user_excluded_in_scaling)
      if self.params.cut_data.d_max:
        reflection.set_flags(reflection['d'] > self.params.cut_data.d_max,
        reflection.flags.user_excluded_in_scaling)
      if self.params.cut_data.partiality_cutoff and 'partiality' in reflection:
        reflection.set_flags(reflection['partiality'] < \
          self.params.cut_data.partiality_cutoff,
          reflection.flags.user_excluded_in_scaling)

    if self.params.cut_data.exclude_image_range:
      if len(self.reflections) == 1:
        start_excl = self.params.cut_data.exclude_image_range[0]
        end_excl = self.params.cut_data.exclude_image_range[1]
        mask1 = start_excl < self.reflections[0]['xyzobs.px.value'].parts()[2]
        mask2 = end_excl > self.reflections[0]['xyzobs.px.value'].parts()[2]
        self.reflections[0].set_flags(mask1 & mask2,
          self.reflections[0].flags.user_excluded_in_scaling)
      else:
        raise Sorry("""exclude_image_range can only be used with one dataset,
        not multiple datasets.""")

  def scale(self):
    """Create the scaling models and perform scaling."""
    self.experiments = create_scaling_model(self.params, self.experiments,
      self.reflections)
    logger.info('\nScaling models have been initialised for all experiments.')
    logger.info('\n' + '='*80 + '\n')

    self.scaler = create_scaler(self.params, self.experiments, self.reflections)
    self.scaler = self.scaling_algorithm(self.scaler)

  def scaled_data_as_miller_array(self, crystal_symmetry, reflection_table=None,
    anomalous_flag=False):
    """Get a scaled miller array from an experiment and reflection table."""
    if not reflection_table:
      joint_table = flex.reflection_table()
      if self.params.scaling_options.only_target or \
        self.params.scaling_options.target_model or \
        self.params.scaling_options.target_mtz:
        self.reflections = self.reflections[:-1]
      for reflection_table in self.reflections:
        #better to just create many miller arrays and join them?
        refl_for_joint_table = flex.reflection_table()
        for col in ['miller_index', 'intensity.scale.value',
          'inverse_scale_factor', 'intensity.scale.variance']:
          refl_for_joint_table[col] = reflection_table[col]
        good_refl_sel = ~reflection_table.get_flags(
          reflection_table.flags.bad_for_scaling, all=False)
        refl_for_joint_table = refl_for_joint_table.select(good_refl_sel)
        joint_table.extend(refl_for_joint_table)
    else:
      good_refl_sel = ~reflection_table.get_flags(
        reflection_table.flags.bad_for_scaling, all=False)
      joint_table = reflection_table.select(good_refl_sel)

    miller_set = miller.set(crystal_symmetry=crystal_symmetry,
      indices=joint_table['miller_index'], anomalous_flag=anomalous_flag)
    i_obs = miller.array(
      miller_set, data=joint_table['intensity.scale.value']/ \
      joint_table['inverse_scale_factor'])
    i_obs.set_observation_type_xray_intensity()
    i_obs.set_sigmas((joint_table['intensity.scale.variance']**0.5)/ \
      joint_table['inverse_scale_factor'])
    i_obs.set_info(
      miller.array_info(source='DIALS', source_type='reflection_tables'))
    return i_obs

  def merging_stats(self):
    """Calculate and print the merging statistics."""
    logger.info('\n'+'='*80+'\n')

    if self.params.output.exclude_on_batch_rmerge:
      self.reflections = exclude_on_batch_rmerge(self.reflections, self.experiments,
        self.params.output.exclude_on_batch_rmerge)

    if self.params.output.exclude_on_image_scale:
      self.reflections = exclude_on_image_scale(self.reflections, self.experiments,
        self.params.output.exclude_on_image_scale)

    self.scaled_miller_array = self.scaled_data_as_miller_array(
      self.experiments[0].crystal.get_crystal_symmetry(), anomalous_flag=False)

    if self.scaled_miller_array.is_unique_set_under_symmetry():
      logger.info(("Dataset doesn't contain any equivalent reflections, \n"
        "no merging statistics can be calculated."))
      return

    make_sub_header("Overall merging statistics (non-anomalous)",
        out=log.info_handle(logger))
    try:
      result = iotbx.merging_statistics.dataset_statistics(
        i_obs=self.scaled_miller_array, n_bins=self.params.output.merging.nbins,
        anomalous=False, sigma_filtering=None, eliminate_sys_absent=False,
        use_internal_variance=self.params.output.use_internal_variance)
      result.show(header=0, out=log.info_handle(logger))
      result.show_estimated_cutoffs(out=log.info_handle(logger))
      self.merging_statistics_result = result
    except RuntimeError:
      raise DialsMergingStatisticsError("Failure during merging statistics calculation")

  def delete_datastructures(self):
    """Delete the data in the scaling datastructures to save RAM before
    combinining datasets for output."""
    del self.scaler
    for experiment in self.experiments:
      for component in experiment.scaling_model.components.iterkeys():
        experiment.scaling_model.components[component] = []
    gc.collect()

  def output(self):
    """Save the experiments json and scaled pickle file."""
    logger.info('\n'+'='*80+'\n')

    if self.params.scaling_options.target_model or \
      self.params.scaling_options.target_mtz or \
      self.params.scaling_options.only_target:
      self.experiments = self.experiments[:-1]
    save_experiments(self.experiments, self.params.output.experiments)

    crystal_symmetry = self.experiments[0].crystal.get_crystal_symmetry() # save for later

    # Now create a joint reflection table. Delete all other data before
    # joining reflection tables - just need experiments for mtz export
    # and a reflection table.
    self.delete_datastructures()

    joint_table = flex.reflection_table()
    for i in range(len(self.reflections)):
      joint_table.extend(self.reflections[i])
      #del reflection_table
      self.reflections[i] = 0
      gc.collect()

    # remove reflections with neg sigma
    sel = joint_table['inverse_scale_factor'] <= 0.0
    n_neg = sel.count(True)
    if n_neg > 0:
      logger.warning(
        'Warning: %s reflections were assigned negative scale factors. \n'
        'It may be best to rerun scaling from this point for an improved model.' % n_neg)
      joint_table.set_flags(sel, joint_table.flags.excluded_for_scaling)

    save_reflections(joint_table, self.params.output.reflections)

    if self.params.output.unmerged_mtz:
      logger.info("\nSaving output to an unmerged mtz file to %s.",
        self.params.output.unmerged_mtz)
      from dials.command_line.export import MTZExporter
      from dials.command_line.export import phil_scope as export_phil_scope
      parser = OptionParser(read_experiments=False, read_reflections=False,
        read_datablocks=False, phil=export_phil_scope)
      params, _ = parser.parse_args(args=[], show_diff_phil=False)
      params.intensity = ['scale']
      params.mtz.partiality_threshold = self.params.cut_data.partiality_cutoff
      params.mtz.hklout = self.params.output.unmerged_mtz
      params.mtz.crystal_name = self.params.output.crystal_name
      if self.params.cut_data.d_min:
        params.mtz.d_min = self.params.cut_data.d_min
      exporter = MTZExporter(params, self.experiments,
        [joint_table])
      exporter.export()

    if self.params.output.merged_mtz:
      logger.info("\nSaving output to a merged mtz file to %s.\n",
        self.params.output.merged_mtz)
      merged_scaled = self.scaled_miller_array.merge_equivalents().array()
      mtz_dataset = merged_scaled.as_mtz_dataset(crystal_name='dials',
        column_root_label='IMEAN') # what does column_root_label do?

      anomalous_scaled = self.scaled_data_as_miller_array(crystal_symmetry,
        reflection_table=joint_table, anomalous_flag=True)
      merged_anom = anomalous_scaled.merge_equivalents(
        use_internal_variance=self.params.output.use_internal_variance).array()
      multiplticies = merged_anom.multiplicities()
      mtz_dataset.add_miller_array(merged_anom, column_root_label='I',
        column_types='KM')
      mtz_dataset.add_miller_array(multiplticies, column_root_label='N',
        column_types='I')
      mtz_file = mtz_dataset.mtz_object()
      mtz_file.set_title('from dials.scale')
      date_str = time.strftime('%d/%m/%Y at %H:%M:%S', time.gmtime())
      mtz_file.add_history('From %s, run on %s' % (dials_version(), date_str))
      mtz_file.write(self.params.output.merged_mtz)

  def scaling_algorithm(self, scaler):
    """The main scaling algorithm."""

    if scaler.id_ == 'target':
      ### FIXME add in quick prescaling round if large scale difference?
      scaler.perform_scaling()

      if scaler.params.scaling_options.only_target or \
        scaler.params.scaling_options.target_model or \
        scaler.params.scaling_options.target_mtz:
        #Do some rounds of targeted scaling and then exit the algorithm.

        if scaler.params.scaling_options.outlier_rejection:
          scaler.outlier_rejection_routine()
          scaler.perform_scaling()

        if scaler.params.scaling_options.full_matrix and (
          scaler.params.scaling_refinery.engine == 'SimpleLBFGS'):
          scaler.perform_scaling(
            engine=scaler.params.scaling_refinery.full_matrix_engine,
            max_iterations=scaler.params.scaling_refinery.full_matrix_max_iterations)

        scaler.expand_scales_to_all_reflections(calc_cov=True)
        if scaler.params.scaling_options.outlier_rejection:
          scaler.round_of_outlier_rejection()
        if scaler.params.weighting.optimise_errors:
          scaler.perform_error_optimisation(update_Ih=False)

        scaler.adjust_variances()

        scaler.clean_reflection_tables()
        return scaler
      # Now pass to a multiscaler ready for next round of scaling.
      scaler.expand_scales_to_all_reflections()
      scaler = MultiScalerFactory.create_from_targetscaler(scaler)

    # From here onwards, scaler should only be a SingleScaler
    # or MultiScaler (not TargetScaler).
    scaler.perform_scaling()

    #Do another round of outlier rejection and then another minimisation.
    rescale = False
    if scaler.params.scaling_options.outlier_rejection:
      scaler.outlier_rejection_routine()
      rescale = True

    if scaler.params.reflection_selection.intensity_choice == 'combine':
      scaler.combine_intensities()
      if scaler.params.scaling_options.outlier_rejection:
        scaler.outlier_rejection_routine()
      rescale = True

    if rescale:
      scaler.perform_scaling()

    # Option to optimise the error model and then do another minimisation.
    if scaler.params.weighting.optimise_errors:
      scaler.error_optimisation_routine()
      scaler.perform_scaling()

    # Now do one round of full matrix minimisation to determine errors.
    if scaler.params.scaling_options.full_matrix and (
      scaler.params.scaling_refinery.engine == 'SimpleLBFGS'):
      scaler.perform_scaling(
        engine=scaler.params.scaling_refinery.full_matrix_engine,
        max_iterations=scaler.params.scaling_refinery.full_matrix_max_iterations)

    # The minimisation has only been done on a subset on the data, so apply the
    # scale factors to the whole reflection table.
    scaler.clear_Ih_table()
    scaler.expand_scales_to_all_reflections(calc_cov=True)
    if scaler.params.scaling_options.outlier_rejection:
      # Note just call the method, not the 'outlier_rejection_routine'
      scaler.round_of_outlier_rejection()

    if scaler.params.weighting.optimise_errors:
      # Note just call the method, not the 'error_optimisation_routine'
      scaler.perform_error_optimisation(update_Ih=False)

    scaler.adjust_variances()
    scaler.clean_reflection_tables()
    return scaler

if __name__ == "__main__":
  try:
    #Parse the command line and flatten reflections, experiments
    usage = '''Usage: dials.scale integrated.pickle integrated_experiments.json
[integrated.pickle(2) integrated_experiments.json(2) ....] [options]'''
    optionparser = OptionParser(usage=usage, read_experiments=True,
      read_reflections=True, read_datablocks=False, phil=phil_scope,
      check_format=False, epilog=help_message)
    params, _ = optionparser.parse_args(show_diff_phil=False)
    if not params.input.experiments or not params.input.reflections:
      optionparser.print_help()
      sys.exit()
    reflections = flatten_reflections(params.input.reflections)
    experiments = flatten_experiments(params.input.experiments)

    if params.cross_validation.cross_validation_mode:
      from dials.algorithms.scaling.cross_validation.cross_validate import \
        cross_validate
      from dials.algorithms.scaling.cross_validation.crossvalidator import \
        DialsScaleCrossValidator

      log.config(verbosity=1, info=params.cross_validation.log,
        debug=params.cross_validation.debug.log)
      logger.info(dials_version())
      diff_phil = optionparser.diff_phil
      if diff_phil.as_str() is not '':
        logger.info('The following parameters have been modified:\n')
        logger.info(diff_phil.as_str())
      diff_phil.objects = [obj for obj in diff_phil.objects if not (
        obj.name == 'input' or obj.name == 'cross_validation')]

      cross_validator = DialsScaleCrossValidator(experiments, reflections)
      cross_validate(params, cross_validator)

      if diff_phil.objects:
        logger.info("\nAdditional configuration for all runs: \n%s", diff_phil.as_str())

      logger.info("Cross validation analysis does not produce scaling output files, rather\n"
        "it gives insight into the dataset. Choose an appropriate parameterisation\n"
        "and rerun scaling without cross_validation_mode.\n")

    else:
      #Set up the log
      log.config(verbosity=1, info=params.output.log,
          debug=params.output.debug.log)
      logger.info(dials_version())
      diff_phil = optionparser.diff_phil.as_str()
      if diff_phil is not '':
        logger.info('The following parameters have been modified:\n')
        logger.info(diff_phil)

      script = Script(params, experiments, reflections)
      script.run()

  except Exception as e:
    halraiser(e)
