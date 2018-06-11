#!/usr/bin/env python
# coding: utf-8
"""
Usage: dials.scale integrated.pickle integrated_experiments.json
[integrated.pickle(2) integrated_experiments.json(2) ....] [options]

This program performs scaling on the input datasets. The default
parameterisation is a physical parameterisation based on that used in the program
Aimless. If multiple input files have been specified, the datasets will be
jointly scaled against a common target of unique reflection intensities.

By default, a scale, decay and absorption correction parameterisation for each
dataset is used. One scaled.pickle and scaled_experiments.json files are output,
which may contain data and scale models from multiple experiments. The
reflection intensities are left unscaled and unmerged in the output, but an
'inverse_scale_factor' and 'inverse_scale_factor_variance' column is added to
the pickle file.

To plot the scale factors determined by this program, one should run:
dials.plot_scaling_models scaled.pickle scaled_experiments.json
"""
from __future__ import absolute_import, division, print_function
import time
import logging
import sys
from copy import deepcopy
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
  calculate_merging_statistics, create_datastructures_for_structural_model
from dials.algorithms.scaling.scaler_factory import create_scaler,\
  MultiScalerFactory
from dials.algorithms.scaling.scaling_utilities import parse_multiple_datasets,\
  select_datasets_on_ids, save_experiments, save_reflections,\
  assign_unique_identifiers


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
    debug_log = dials.scale.debug.log
      .type = str
      .help = "The debug log filename"
    calculate_individual_merging_stats = False
      .type = bool
      .help = "Option to calculate merging stats for the individual datasets."
    plot_merging_stats = False
      .type = bool
      .help = "Option to switch on plotting of merging stats."
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
    use_internal_variance = False
      .type = bool
      .help = "Option to use internal spread of the intensities when merging
              reflection groups and calculating sigI, rather than using the
              sigmas of the individual reflections."
  }
  include scope dials.algorithms.scaling.scaling_options.phil_scope
  include scope dials.algorithms.scaling.scaling_refiner.scaling_refinery_phil_scope
''', process_includes=True)


class Script(object):
  """Main script to run the scaling algorithm."""

  def __init__(self, params, experiments, reflections):
    self.params = params
    self.experiments = experiments
    self.reflections = reflections
    self.scaler = None
    self.minimised = None
    self.scaled_miller_array = None
    self.dataset_ids = []

  def run(self, save_data=True):
    """Run the scaling script."""
    start_time = time.time()
    self.prepare_input()
    self.scale()
    self.merging_stats()
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
      if experiment.crystal.get_space_group() != s_g_1:
        raise Sorry('experiments have different space groups and cannot be'
          'scaled together, please reanalyse the data so that the space groups'
          'are consistent or manually specify a space group.'
          'Alternatively, some datasets can be excluded using exclude_datasets=')

    if self.params.scaling_options.target_model:
      logger.info("Extracting data from structural model.")
      exp, reflections = create_datastructures_for_structural_model(
        self.reflections, self.experiments, self.params.scaling_options.target_model)
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
    self.minimised = self.scaling_algorithm(self.scaler)

  def scaled_data_as_miller_array(self, experiment, reflection_table,
    anomalous_flag=False):
    """Get a scaled miller array from an experiment and reflection table."""
    bad_refl_sel = reflection_table.get_flags(
      reflection_table.flags.bad_for_scaling, all=False)
    r_t = reflection_table.select(~bad_refl_sel)
    miller_set = miller.set(
      crystal_symmetry=experiment.crystal.get_crystal_symmetry(),
      indices=r_t['miller_index'], anomalous_flag=anomalous_flag)
    i_obs = miller.array(
      miller_set, data=r_t['intensity']/r_t['inverse_scale_factor'])
    i_obs.set_observation_type_xray_intensity()
    i_obs.set_sigmas((r_t['variance']**0.5)/r_t['inverse_scale_factor'])
    i_obs.set_info(
      miller.array_info(source='DIALS', source_type='reflection_tables'))
    return i_obs


  def merging_stats(self):
    """Calculate and print the merging statistics."""
    logger.info('\n'+'='*80+'\n')
    # Calculate merging stats.
    plot_labels = []


    if self.params.output.calculate_individual_merging_stats and (
      len(set(self.minimised.reflection_table['id'])) > 1):
      results, scaled_ids = calculate_merging_statistics(
        self.minimised.reflection_table, self.experiments,
        use_internal_variance=self.params.output.use_internal_variance)
      for result, data_id in zip(results, scaled_ids):
        make_sub_header("Merging statistics for dataset " + str(data_id),
          out=log.info_handle(logger))
        result.show(header=0, out=log.info_handle(logger))
        result.show_estimated_cutoffs(out=log.info_handle(logger))
        plot_labels.append('Dataset ' + str(data_id))

    self.scaled_miller_array = self.scaled_data_as_miller_array(self.experiments[0],
      self.minimised.reflection_table, anomalous_flag=False)

    make_sub_header("Overall merging statistics (non-anomalous)",
        out=log.info_handle(logger))
    result = iotbx.merging_statistics.dataset_statistics(
      i_obs=self.scaled_miller_array, n_bins=20, anomalous=False,
      sigma_filtering=None, eliminate_sys_absent=False,
      use_internal_variance=self.params.output.use_internal_variance)
    result.show(header=0, out=log.info_handle(logger))
    result.show_estimated_cutoffs(out=log.info_handle(logger))
    plot_labels.append('Overall dataset')

    # Plot merging stats if requested.
    if self.params.output.plot_merging_stats:
      from xia2.command_line.compare_merging_stats import plot_merging_stats
      plot_merging_stats(results, labels=plot_labels)

  def output(self):
    """Save the experiments json and scaled pickle file."""
    logger.info('\n'+'='*80+'\n')

    if self.params.scaling_options.target_model:
      self.experiments = self.experiments[:-1]

    self.minimised.reflection_table['intensity.scale.value'] = deepcopy(
      self.minimised.reflection_table['intensity'])
    self.minimised.reflection_table['intensity.scale.variance'] = deepcopy(
      self.minimised.reflection_table['variance'])

    if self.params.output.unmerged_mtz:
      logger.info("\nSaving output to an unmerged mtz file to %s.",
        self.params.output.unmerged_mtz)
      from dials.command_line.export import MTZExporter
      from dials.command_line.export import phil_scope as export_phil_scope
      parser = OptionParser(read_experiments=False, read_reflections=False,
        read_datablocks=False, phil=export_phil_scope)
      params, _ = parser.parse_args(args=[], show_diff_phil=False)
      params.mtz.apply_scales = True
      params.mtz.hklout = self.params.output.unmerged_mtz
      if self.params.scaling_options.integration_method == 'sum':
        params.mtz.ignore_profile_fitting = True #to make it export summation
      exporter = MTZExporter(params, self.experiments,
        [self.minimised.reflection_table])
      exporter.export()

    if self.params.output.merged_mtz:
      logger.info("\nSaving output to a merged mtz file to %s.\n",
        self.params.output.merged_mtz)
      merged_scaled = self.scaled_miller_array.merge_equivalents().array()
      mtz_dataset = merged_scaled.as_mtz_dataset(crystal_name='dials',
        column_root_label='IMEAN') # what does column_root_label do?
      anomalous_scaled = self.scaled_data_as_miller_array(self.experiments[0],
        self.minimised.reflection_table, anomalous_flag=True)
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

    save_experiments(self.experiments, self.params.output.experiments)
    save_reflections(self.minimised, self.params.output.reflections)

    '''if params.output.plot_scaling_models:
      from dials.command_line.plot_scaling_models import plot_scaling_models
      plot_scaling_models(params.output.scaled, params.output.experiments)'''

  def scaling_algorithm(self, scaler):
    """The main scaling algorithm."""

    if scaler.id_ == 'target':
      ### FIXME add in quick prescaling round if large scale difference?
      scaler.perform_scaling()

      if scaler.params.scaling_options.only_target or (
        scaler.params.scaling_options.target_model):
        #Do some rounds of targeted scaling and then exit the algorithm.
        #scaler.expand_scales_to_all_reflections()
        # Do another round so that more suitable weights are used.
        #scaler.select_reflections_for_scaling()
        #scaler.perform_scaling() ##Note - were these helping?

        if scaler.params.scaling_options.full_matrix and (
          scaler.params.scaling_refinery.engine == 'SimpleLBFGS'):
          scaler.perform_scaling(
            engine=scaler.params.scaling_refinery.full_matrix_engine,
            max_iterations=scaler.params.scaling_refinery.full_matrix_max_iterations)

        scaler.expand_scales_to_all_reflections(calc_cov=True)

        if scaler.params.weighting.optimise_errors:
          scaler.perform_error_optimisation()

        scaler.adjust_variances()

        scaler.join_multiple_datasets()
        return scaler
      # Now pass to a multiscaler ready for next round of scaling.
      scaler.expand_scales_to_all_reflections()
      scaler = MultiScalerFactory.create_from_targetscaler(scaler)

    # From here onwards, scaler should only be a SingleScaler
    # or MultiScaler (not TargetScaler).
    scaler.perform_scaling()

    #Do another round of outlier rejection and then another minimisation.
    if scaler.params.scaling_options.outlier_rejection:
      scaler.outlier_rejection_routine()
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
    scaler.expand_scales_to_all_reflections(calc_cov=True)
    if scaler.params.scaling_options.outlier_rejection:
      # Note just call the method, not the 'outlier_rejection_routine'
      scaler.round_of_outlier_rejection()

    if scaler.params.weighting.optimise_errors:
      # Note just call the method, not the 'error_optimisation_routine'
      scaler.perform_error_optimisation()

    scaler.adjust_variances()

    if scaler.id_ == 'multi':
      scaler.join_multiple_datasets()
    return scaler

if __name__ == "__main__":
  try:
    #Parse the command line and flatten reflections, experiments
    optionparser = OptionParser(usage=__doc__.strip(), read_experiments=True,
      read_reflections=True, read_datablocks=False, phil=phil_scope,
      check_format=False)
    params, _ = optionparser.parse_args(show_diff_phil=False)
    if not params.input.experiments or not params.input.reflections:
      optionparser.print_help()
      sys.exit()
    reflections = flatten_reflections(params.input.reflections)
    experiments = flatten_experiments(params.input.experiments)

    #Set up the log
    log.config(verbosity=1, info=params.output.log,
        debug=params.output.debug_log)
    logger.info(dials_version())
    diff_phil = optionparser.diff_phil.as_str()
    if diff_phil is not '':
      logger.info('The following parameters have been modified:\n')
      logger.info(diff_phil)

    script = Script(params, experiments, reflections)
    script.run()

  except Exception as e:
    halraiser(e)
