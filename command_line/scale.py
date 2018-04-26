#!/usr/bin/env python
# coding: utf-8
"""
Usage: dials_scratch.scale integrated.pickle integrated_experiments.json
[integrated.pickle(2) integrated_experiments.json(2) ....] [options]

This program performs scaling on the input datasets. The default
parameterisationis a physical parameterisation based on that used in the program
Aimless. If multiple input files have been specified, the datasets will be
jointly scaled against a common target of unique reflection intensities.

By default, a scale, decay and absorption correction parameterisation for each
dataset is used. One scaled.pickle and scaled_experiments.json files are output,
which may contain data and scale models from multiple experiments. The
reflection intensities are left unscaled and unmerged in the output, but an
'inverse_scale_factor' and 'inverse_scale_factor_variance' column is added to
the pickle file.

To plot the scale factors determined by this program, one should run:
dials_scratch.plot_scaling_models scaled.pickle scaled_experiments.json
"""
from __future__ import absolute_import, division, print_function
import time
import logging
#import libtbx.load_env
from libtbx import phil
from libtbx.utils import Sorry
from libtbx.str_utils import make_sub_header
from cctbx.miller import crystal
from dials.util import halraiser, log
from dials.array_family import flex
from dials.util.options import OptionParser, flatten_reflections,\
  flatten_experiments
from dials.algorithms.scaling.scaling_library import create_scaling_model,\
  calculate_merging_statistics, calculate_single_merging_stats
from dials.algorithms.scaling.scaler_factory import create_scaler,\
  MultiScalerFactory
from dials.algorithms.scaling.scaling_utilities import (
  parse_multiple_datasets, save_experiments, save_reflections)


logger = logging.getLogger('dials')

phil_scope = phil.parse('''
  debug = False
    .type = bool
    .help = "Output additional debugging information"
  model = *physical array KB
      .type = choice
      .help = "Set scaling model to be applied to input datasets without
               an existing model. "
  output {
    log = dials_scratch.scaling.log
      .type = str
      .help = "The log filename"
    debug_log = dials_scratch.scaling.debug.log
      .type = str
      .help = "The debug log filename"
    plot_merging_stats = False
      .type = bool
      .help = "Option to switch on plotting of merging stats."
    plot_scaling_models = False
      .type = bool
      .help = "Option to switch on plotting of the scaling models determined."
    experiments = "scaled_experiments.json"
      .type = str
      .help = "Option to set filepath for output json."
    scaled = "scaled.pickle"
      .type = str
      .help = "Option to set filepath for output pickle file of scaled
               intensities."
  }
  include scope dials.algorithms.scaling.scaling_options.phil_scope
  include scope dials.algorithms.scaling.scaling_refiner.scaling_refinery_phil_scope
''', process_includes=True)


class Script(object):
  """Main script to run the scaling algorithm."""

  def __init__(self):
    optionparser = OptionParser(usage=__doc__.strip(), read_experiments=True,
      read_reflections=True, read_datablocks=False, phil=phil_scope,
      check_format=False)
    self.params, _ = optionparser.parse_args(show_diff_phil=False)
    self.scaler = None
    self.minimised = None

    log.config(verbosity=1, info=self.params.output.log,
      debug=self.params.output.debug_log)

    if not self.params.input.experiments or not self.params.input.reflections:
      optionparser.print_help()
      return

    from dials.util.version import dials_version
    logger.info(dials_version())

    diff_phil = optionparser.diff_phil.as_str()
    if diff_phil is not '':
      logger.info('The following parameters have been modified:\n')
      logger.info(diff_phil)

    self.reflections = flatten_reflections(self.params.input.reflections)
    self.experiments = flatten_experiments(self.params.input.experiments)

    if self.params.scaling_options.target_intensities:
      target_params, _ = optionparser.parse_args(
        [self.params.scaling_options.target_intensities], show_diff_phil=False)
      target_reflections = flatten_reflections(target_params.input.reflections)
      self.reflections.append(target_reflections[0])

    if self.params.scaling_options.target_intensities:
      logger.info('*'*36 + ' WARNING ' + '*'*36 + '\n'
      'For scaling against target calculated dataset, it is assumed \n'
      'that this is the last pickle file passed in and that no corresponding \n'
      'experiments file is passed in.\n')
      self.experiments.append(self.experiments[0]) # Assume correct space group for now.
      self.reflections[-1]['id'] = flex.int(self.reflections[-1].size(),
        len(self.reflections))

    if len(self.experiments) != 1:
      logger.info('Checking for the existence of a reflection table \n'
        'containing multiple scaled datasets \n')
      self.reflections = parse_multiple_datasets(self.reflections)
      logger.info("Found %s reflection tables in total.", len(self.reflections))
      logger.info("Found %s experiments in total.", len(self.experiments))
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
            'are consistent or manually specify a space group.')

    if len(self.experiments) != len(self.reflections):
      raise Sorry("Mismatched number of experiments and reflection tables found.")

  def run(self):
    """Run the scaling script."""
    start_time = time.time()
    self.prepare_input()
    self.scale()
    self.merging_stats()
    self.output()
    # All done!
    finish_time = time.time()
    logger.info("\nTotal time taken: {0:.4f}s ".format(finish_time - start_time))
    logger.info('\n'+'='*80+'\n')


  def prepare_input(self):
    """Perform any cutting of the dataset prior to creating scaling models."""
    for reflection in self.reflections:
      reflection.set_flags(flex.bool(reflection.size(), False),
        reflection.flags.user_excluded_in_scaling)
      if self.params.cut_data.max_resolution:
        reflection.set_flags(reflection['d'] < self.params.cut_data.max_resolution,
        reflection.flags.user_excluded_in_scaling)
      if self.params.cut_data.min_resolution:
        reflection.set_flags(reflection['d'] > self.params.cut_data.min_resolution,
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


  def merging_stats(self):
    """Calculate and print the merging statistics."""
    logger.info('\n'+'='*80+'\n')
    # Calculate merging stats.
    plot_labels = []
    results, scaled_ids = calculate_merging_statistics(
      self.minimised.reflection_table, self.experiments)
    if len(results) == 1:
      logger.info("")
      results[0].show()#out=log.info_handle(logger))
      results[0].show_estimated_cutoffs()#out=log.info_handle(logger))
      plot_labels.append('Single dataset ')
    else:
      for result, data_id in zip(results, scaled_ids):
        make_sub_header("Merging statistics for dataset " + str(data_id),
          out=log.info_handle(logger))
        result.show(header=0)#out=log.info_handle(logger))
        result.show_estimated_cutoffs()
        plot_labels.append('Dataset ' + str(data_id))
      make_sub_header("Merging statistics for combined datasets",
        out=log.info_handle(logger))
      result = calculate_single_merging_stats(
        self.minimised.reflection_table, self.experiments[0])
      result.show(header=0)#out=log.info_handle(logger))
      result.show_estimated_cutoffs()
      plot_labels.append('Combined datasets')

    # Plot merging stats if requested.
    if self.params.output.plot_merging_stats:
      from xia2.command_line.compare_merging_stats import plot_merging_stats
      plot_merging_stats(results, labels=plot_labels)

  def output(self):
    """Save the experiments json and scaled pickle file."""
    logger.info('\n'+'='*80+'\n')
    if self.params.scaling_options.target_intensities:
      self.experiments = self.experiments[:-1]
    save_experiments(self.experiments, self.params.output.experiments)
    save_reflections(self.minimised, self.params.output.scaled)

    '''if params.output.plot_scaling_models:
      from dials_scratch.command_line.plot_scaling_models import plot_scaling_models
      plot_scaling_models(params.output.scaled, params.output.experiments)'''

  def scaling_algorithm(self, scaler):
    """The main scaling algorithm."""

    if scaler.id_ == 'target':
      scaler.perform_scaling()

      # The minimisation has only been done on a subset on the data, so apply the
      # scale factors to the whole reflection table.
      if scaler.params.scaling_options.only_target or (
        scaler.params.scaling_options.target_intensities):
        scaler.expand_scales_to_all_reflections()

        # Do another round so that more suitable weights are used.
        scaler.select_reflections_for_scaling()
        scaler.perform_scaling()

        #Need to add in a full matrix round to allow calculation of variances.

        scaler.expand_scales_to_all_reflections()

        if scaler.params.weighting.optimise_error_model:
          scaler.perform_error_optimisation()
          '''scaler.apply_error_model_to_variances()'''

        scaler.join_multiple_datasets()
        return scaler
      # Now pass to a multiscaler ready for next round of scaling.
      scaler.expand_scales_to_all_reflections()
      scaler = MultiScalerFactory.create_from_targetscaler(scaler)

    # From here onwards, scaler should only be a SingleScaler
    # or MultiScaler (not TargetScaler).
    scaler.perform_scaling()

    # Option to optimise the error model and then do another minimisation.
    if scaler.params.weighting.optimise_error_model:
      scaler.expand_scales_to_all_reflections()
      scaler.perform_error_optimisation()
      scaler.select_reflections_for_scaling()
      scaler.perform_scaling()

    # Now do one round of full matrix minimisation to determine errors.
    if scaler.params.scaling_options.full_matrix_round and (
      scaler.params.scaling_refinery.engine == 'SimpleLBFGS'):
      scaler.perform_scaling(
        engine=scaler.params.scaling_refinery.full_matrix_engine,
        max_iterations=scaler.params.scaling_refinery.full_matrix_max_iterations)

    # The minimisation has only been done on a subset on the data, so apply the
    # scale factors to the whole reflection table.
    scaler.expand_scales_to_all_reflections(calc_cov=True)

    if scaler.params.weighting.optimise_error_model:
      scaler.perform_error_optimisation()
      #scaler.apply_error_model_to_variances()

    if scaler.id_ == 'multi':
      scaler.join_multiple_datasets()
    return scaler

if __name__ == "__main__":
  try:
    script = Script()
    script.run()
  except Exception as e:
    halraiser(e)
