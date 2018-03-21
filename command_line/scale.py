#!/usr/bin/env python
# coding: utf-8
"""
Usage: dials_scratch.scale integrated.pickle integrated_experiments.json
[integrated.pickle(2) integrated_experiments.json(2) ....] [options]

This program performs scaling on the input datasets. The default
parameterisationis a physical parameterisation based on that used in the program
Aimless. Ifmultiple input files have been specified, the datasets will be
jointly scaledagainst a common target of unique reflection intensities.

By default, a scale, decay and absorption correction parameterisation for each
dataset is used. One scaled.pickle and scaled_experiments.json files are output,
which may contain data and scale models from multiple experiments. The
reflectionintensities are left unscaled and unmerged in the output, but an
'inverse_scale_factor' and 'inverse_scale_factor_variance' column is added to
the pickle file.

To plot the scale factors determined by this program, one should run:
dials_scratch.plot_scaling_models scaled.pickle scaled_experiments.json
"""
from __future__ import absolute_import, division, print_function
import time
import logging
import sys
import libtbx.load_env
from libtbx import phil
from libtbx.utils import Sorry
from libtbx.str_utils import make_sub_header
from dials.util import halraiser
from dials.util.options import OptionParser, flatten_reflections,\
  flatten_experiments
from dials.algorithms.scaling.scaling_refiner import scaling_refinery,\
  error_model_refinery
from dials.algorithms.scaling.model.scaling_model_factory import \
  create_scaling_model
from dials.algorithms.scaling.scaler_factory import create_scaler,\
  MultiScalerFactory
from dials.algorithms.scaling import scaler as scaler_module
from dials.algorithms.scaling.parameter_handler import create_apm
from dials.algorithms.scaling.target_function import ScalingTarget,\
  ScalingTargetFixedIH
from dials.algorithms.scaling.error_model.error_model import \
  BasicErrorModel
from dials.algorithms.scaling.error_model.error_model_target import \
  ErrorModelTarget
from dials.algorithms.scaling.scaling_utilities import (
  parse_multiple_datasets, save_experiments, save_reflections)


start_time = time.time()
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

def main(argv):
  """Main script to run the scaling algorithm."""

  optionparser = OptionParser(usage=__doc__.strip(), read_experiments=True,
    read_reflections=True, read_datablocks=False, phil=phil_scope,
    check_format=False)
  params, _ = optionparser.parse_args(argv, show_diff_phil=False)

  from dials.util import log
  log.config(verbosity=1, info=params.output.log, debug=params.output.debug_log)

  if not params.input.experiments or not params.input.reflections:
    optionparser.print_help()
    return

  from dials.util.version import dials_version
  logger.info(dials_version())

  diff_phil = optionparser.diff_phil.as_str()
  if diff_phil is not '':
    logger.info('The following parameters have been modified:\n')
    logger.info(diff_phil)

  reflections = flatten_reflections(params.input.reflections)
  experiments = flatten_experiments(params.input.experiments)

  if params.scaling_options.target_intensities:
    from dials.array_family import flex
    target_params, _ = optionparser.parse_args(
      [params.scaling_options.target_intensities], show_diff_phil=False)
    target_reflections = flatten_reflections(target_params.input.reflections)
    reflections.append(target_reflections[0])

  # If scaling against a target calculated dataset, it is assumed that this
  # is the last pickle file passed in and that no corresponding experiments file
  # is passed in.
  if params.scaling_options.target_intensities:
    logger.info('*'*36 + ' WARNING ' + '*'*36 + '\n'
    'For scaling against target calculated dataset, it is assumed \n'
    'that this is the last pickle file passed in and that no corresponding \n'
    'experiments file is passed in.\n')
    experiments.append(experiments[0]) # Assume correct space group.
    reflections[-1]['id'] = flex.int([len(reflections)] * reflections[-1].size())

  if len(experiments) != 1:
    logger.info(('Checking for the existence of a reflection table {sep}'
      'containing multiple scaled datasets {sep}').format(sep='\n'))
    reflections = parse_multiple_datasets(reflections)
    logger.info("Found %s reflection tables in total." % len(reflections))
    logger.info("Found %s experiments in total." % len(experiments))
    s_g_1 = experiments[0].crystal.get_space_group()
    for experiment in experiments:
      if experiment.crystal.get_space_group() != s_g_1:
        raise Sorry("""experiments have different space groups and cannot be
        scaled together, please reanalyse the data so that the space groups
        are consistent""")

  if len(experiments) != len(reflections):
     raise Sorry("Mismatched number of experiments and reflection tables found.")

  # Perform any cutting of the dataset before creating scaling models.
  for reflection in reflections:
    reflection.set_flags(flex.bool([False]*reflection.size()),
      reflection.flags.user_excluded_in_scaling)
    if params.cut_data.max_resolution:
      reflection.set_flags(reflection['d'] < params.cut_data.max_resolution,
      reflection.flags.user_excluded_in_scaling)
    if params.cut_data.min_resolution:
      reflection.set_flags(reflection['d'] > params.cut_data.min_resolution,
      reflection.flags.user_excluded_in_scaling)
  if params.cut_data.exclude_image_range:
    if len(reflections) == 1:
      start_excl = params.cut_data.exclude_image_range[0]
      end_excl = params.cut_data.exclude_image_range[1]
      mask1 = start_excl < reflections[0]['xyzobs.px.value'].parts()[2]
      mask2 = end_excl > reflections[0]['xyzobs.px.value'].parts()[2]
      reflections[0].set_flags(mask1 & mask2,
        reflections[0].flags.user_excluded_in_scaling)
    else:
      raise Sorry("""exclude_image_range can only be used with one dataset,
      not multiple datasets.""")
  

  # First create the scaling model if it didn't already exist in the
  # experiments files.
  experiments = create_scaling_model(params, experiments, reflections)
  logger.info('\nScaling models have been initialised for all experiments.')
  logger.info('\n' + '='*80 + '\n')

  # Now create the scaler and do the scaling.
  scaler = create_scaler(params, experiments, reflections)
  minimised = scaling_algorithm(scaler)

  for experiment in experiments:
    experiment.scaling_model.set_scaling_model_as_scaled()

  logger.info('\n'+'='*80+'\n')
  # Calculate merging stats.
  results, scaled_ids = minimised.calc_merging_statistics()
  plot_labels = []
  for i, result in enumerate(results):
    if len(results) == 1:
      logger.info("")
      result.show()#out=log.info_handle(logger))
      plot_labels.append('Single dataset ')
    else:
      if i < len(results) - 1:
        make_sub_header("Merging statistics for dataset " + str(scaled_ids[i]),
          out=log.info_handle(logger))
        result.show(header=0)#out=log.info_handle(logger))
        plot_labels.append('Dataset ' + str(scaled_ids[i]))
      else:
        make_sub_header("Merging statistics for combined datasets",
          out=log.info_handle(logger))
        result.show(header=0)#out=log.info_handle(logger))
        plot_labels.append('Combined datasets')

  # Plot merging stats if requested.
  if params.output.plot_merging_stats:
    from xia2.command_line.compare_merging_stats import plot_merging_stats
    plot_merging_stats(results, labels=plot_labels)
  
  logger.info('\n'+'='*80+'\n')
  # Save scaled_experiments.json and scaled.pickle files.
  if params.scaling_options.target_intensities:
    experiments = experiments[:-1]
  save_experiments(experiments, params.output.experiments)
  minimised.clean_reflection_table()
  save_reflections(minimised.reflection_table, params.output.scaled)

  '''if params.output.plot_scaling_models:
    from dials_scratch.command_line.plot_scaling_models import plot_scaling_models
    plot_scaling_models(params.output.scaled, params.output.experiments)'''

  # All done!
  finish_time = time.time()
  logger.info("\nTotal time taken: {0:.4f}s ".format(finish_time - start_time))
  logger.info('\n'+'='*80+'\n')

def perform_scaling(scaler, target_type=ScalingTarget):
  """Perform a complete minimisation based on the current state."""
  apm_factory = create_apm(scaler)
  for _ in range(apm_factory.n_cycles):
    apm = apm_factory.make_next_apm()
    refinery = scaling_refinery(engine=scaler.params.scaling_refinery.engine,
      target=target_type(scaler, apm), prediction_parameterisation=apm,
      max_iterations=scaler.params.scaling_refinery.max_iterations)
    refinery.run()
    scaler = refinery.return_scaler()
    logger.info('\n'+'='*80+'\n')
  return scaler

def perform_error_optimisation(scaler):
  """Perform an optimisation of the sigma values."""
  if isinstance(scaler, scaler_module.MultiScalerBase):
    for s_scaler in scaler.single_scalers:
      refinery = error_model_refinery(engine='SimpleLBFGS',
        target=ErrorModelTarget(BasicErrorModel(s_scaler.Ih_table)),
        max_iterations=100)
      refinery.run()
      error_model = refinery.return_error_manager()
      s_scaler.update_error_model(error_model.refined_parameters)
    scaler.Ih_table.update_weights_from_error_models()
  else:
    refinery = error_model_refinery(engine='SimpleLBFGS',
      target=ErrorModelTarget(BasicErrorModel(scaler.Ih_table)),
      max_iterations=100)
    refinery.run()
    error_model = refinery.return_error_manager()
    scaler.update_error_model(error_model.refined_parameters)
  return scaler

def scaling_algorithm(scaler):
  """The main scaling algorithm."""

  if scaler.id_ == 'target':
    scaler = perform_scaling(scaler, target_type=ScalingTargetFixedIH)

    # The minimisation has only been done on a subset on the data, so apply the
    # scale factors to the whole reflection table.
    if scaler.params.scaling_options.only_target or (
      scaler.params.scaling_options.target_intensities):
      scaler.expand_scales_to_all_reflections()

      # Do another round so that more suitable weights are used.
      scaler.select_reflections_for_scaling()
      scaler = perform_scaling(scaler, target_type=ScalingTargetFixedIH)
      scaler.expand_scales_to_all_reflections()

      if scaler.params.weighting.optimise_error_model:
        scaler = perform_error_optimisation(scaler)
        '''scaler.apply_error_model_to_variances()'''

      scaler.join_multiple_datasets()
      return scaler
    # Now pass to a multiscaler ready for next round of scaling.
    scaler.expand_scales_to_all_reflections(calc_cov=False)
    scaler = MultiScalerFactory.create_from_targetscaler(scaler)

  # From here onwards, scaler should only be a SingleScaler
  # or MultiScaler (not TargetScaler).
  scaler = perform_scaling(scaler)

  # Option to optimise the error model and then do another minimisation.
  if scaler.params.weighting.optimise_error_model:
    scaler.expand_scales_to_all_reflections()
    scaler = perform_error_optimisation(scaler)
    scaler.select_reflections_for_scaling()
    scaler = perform_scaling(scaler)

  # Now do one round of full matrix minimisation to determine errors.
  if scaler.params.scaling_options.full_matrix_round:
    apm_factory = create_apm(scaler)
    for _ in range(apm_factory.n_cycles):
      apm = apm_factory.make_next_apm()
      refinery = scaling_refinery(
        engine=scaler.params.scaling_refinery.full_matrix_engine,
        target=ScalingTarget(scaler, apm), prediction_parameterisation=apm,
        max_iterations=scaler.params.scaling_refinery.full_matrix_max_iterations)
      try:
        refinery.run()
      except:
        raise Sorry(""" Unable to do a round of full matrix minimisation needed
        to determine scale factor uncertainties. This is likely due to
        overparameterisation, and can often be avoided by reducing the number
        of model parameters through the command line options, or by removing
        model components (e.g. absorption_term=False). As a last resort,
        the full matrix round can be turned off with full_matrix_round=False.""")
      scaler = refinery.return_scaler()
      logger.info('\n'+'='*80+'\n')

  # The minimisation has only been done on a subset on the data, so apply the
  # scale factors to the whole reflection table.
  scaler.expand_scales_to_all_reflections()

  if scaler.params.weighting.optimise_error_model:
    scaler = perform_error_optimisation(scaler)
    '''scaler.apply_error_model_to_variances()'''

  if isinstance(scaler, scaler_module.MultiScalerBase):
    scaler.join_multiple_datasets()
  return scaler

if __name__ == "__main__":
  try:
    sys.exit(main(sys.argv[1:]))
  except Exception as e:
    halraiser(e)
