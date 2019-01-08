'''
Collection of factories for creating the scalers.
'''
import logging
from libtbx.utils import Sorry
from dials.array_family import flex
from dials.algorithms.scaling.scaler import MultiScaler, TargetScaler,\
  SingleScaler
from dials.algorithms.scaling.scaling_utilities import quasi_normalisation, \
  Reasons, BadDatasetForScalingException
from dials.algorithms.scaling.scaling_library import choose_scaling_intensities
from dials.algorithms.scaling.outlier_rejection import reject_outliers
logger = logging.getLogger('dials')

def create_scaler(params, experiments, reflections):
  """Read an experimentlist and list of reflection tables and return
    an appropriate scaler. Requires experiment identifiers are correctly set in
    the experiments and reflections."""
  if len(reflections) == 1:
    scaler = SingleScalerFactory.create(params, experiments[0], reflections[0])
  else:
    is_scaled_list = is_scaled(experiments)
    n_scaled = is_scaled_list.count(True)
    if (params.scaling_options.target_cycle and n_scaled < len(reflections)
        and n_scaled > 0): #if only some scaled and want to do targeted scaling
      scaler = TargetScalerFactory.create(params, experiments, reflections,
        is_scaled_list)
    elif len(reflections) > 1: #else just make one multiscaler for all refls
      scaler = MultiScalerFactory.create(params, experiments, reflections)
    else:
      raise Sorry("no reflection tables found to create a scaler")
  return scaler

def is_scaled(experiments):
  'helper function to return a boolean list of whether experiments are scaled'
  is_already_scaled = []
  for experiment in experiments:
    if experiment.scaling_model.is_scaled:
      is_already_scaled.append(True)
    else:
      is_already_scaled.append(False)
  return is_already_scaled

class ScalerFactory(object):
  """Base class for Scaler Factories"""
  @classmethod
  def filter_bad_reflections(cls, reflections):
    """Initial filter to select integrated reflections."""
    reasons = Reasons()
    mask = ~reflections.get_flags(reflections.flags.integrated, all=False)
    reasons.add_reason('not integrated', mask.count(True))
    if 'd' in reflections:
      d_mask = reflections['d'] <= 0.0
      reasons.add_reason('bad d-value', d_mask.count(True))
      mask = mask | d_mask
    reflections.set_flags(mask, reflections.flags.excluded_for_scaling)
    return reflections, reasons

  @classmethod
  def ensure_experiment_identifier(cls, params, experiment, reflection_table):
    """Check for consistent experiment identifier, and if not set then set it
    using scaled_id."""
    id_vals = reflection_table.experiment_identifiers().values()
    assert experiment.identifier in id_vals, (experiment.identifier, list(id_vals))
    assert len(id_vals) == 1, list(id_vals)
    if params.scaling_options.verbosity > 1:
      logger.info('The experiment identifier for this dataset is %s',
        experiment.identifier)

class SingleScalerFactory(ScalerFactory):
  """Factory for creating a scaler for a single dataset."""

  @classmethod
  def create(cls, params, experiment, reflection_table, for_multi=False):
    """Perform reflection_table preprocessing and create a SingleScaler."""

    cls.ensure_experiment_identifier(params, experiment, reflection_table)

    if params.scaling_options.verbosity > 1:
      logger.info('Preprocessing data for scaling. The id assigned to this \n'
        'dataset is %s, and the scaling model type being applied is %s. \n' %
        (reflection_table.experiment_identifiers().values()[0], experiment.scaling_model.id_))

    reflection_table, reasons = cls.filter_bad_reflections(reflection_table)
    excluded_for_scaling = reflection_table.get_flags(
      reflection_table.flags.excluded_for_scaling)
    user_excluded = reflection_table.get_flags(
      reflection_table.flags.user_excluded_in_scaling)
    n_excluded = (excluded_for_scaling | user_excluded).count(True)
    if n_excluded == reflection_table.size():
      logger.info("All reflections were determined to be unsuitable for scaling.")
      logger.info(reasons)
      raise BadDatasetForScalingException("""Unable to use this dataset for scaling""")
    elif params.scaling_options.verbosity > 1:
      logger.info('%s reflections not suitable for scaling', n_excluded)
      logger.info(reasons)
    if not 'inverse_scale_factor' in reflection_table:
      reflection_table['inverse_scale_factor'] = flex.double(
        reflection_table.size(), 1.0)

    reflection_table = choose_scaling_intensities(reflection_table,
      params.reflection_selection.intensity_choice)

    reflection_table = cls.filter_outliers(reflection_table, experiment,
      params)

    return SingleScaler(params, experiment, reflection_table, for_multi)

  @classmethod
  def filter_outliers(cls, reflections, experiment, params):
    """Calculate normalised E2 values and perform outlier rejection."""
    reflections = quasi_normalisation(reflections, experiment)
    '''if params.scaling_options.outlier_rejection:
      reflections = reject_outliers([reflections],
        experiment.crystal.get_space_group(),
        params.scaling_options.outlier_rejection,
        params.scaling_options.outlier_zmax)[0]'''
    return reflections

class NullScalerFactory(ScalerFactory):
  'Factory for creating null scaler'
  @classmethod
  def create(cls, params, experiment, reflection_table):
    """Return Null Scaler."""
    from dials.algorithms.scaling.scaler import NullScaler
    logger.info('Preprocessing target dataset for scaling. \n')
    reflection_table, reasons = cls.filter_bad_reflections(reflection_table)
    variance_mask = reflection_table['variance'] <= 0.0
    reflection_table.set_flags(variance_mask,
      reflection_table.flags.excluded_for_scaling)
    logger.info('%s reflections not suitable for scaling\n', reflection_table.get_flags(
        reflection_table.flags.excluded_for_scaling).count(True))
    logger.info(reasons)
    cls.ensure_experiment_identifier(params, experiment, reflection_table)
    return NullScaler(params, experiment, reflection_table)

class MultiScalerFactory(object):
  'Factory for creating a scaler for multiple datasets'
  @classmethod
  def create(cls, params, experiments, reflections):
    '''create a list of single scalers to pass to a MultiScaler.'''
    single_scalers = []
    offset = 0
    for i in range(len(reflections)):
      # Remove bad datasets that literally have no integrated reflections
      try:
        scaler = SingleScalerFactory.create(
          params, experiments[i-offset], reflections[i-offset], for_multi=True)
        single_scalers.append(scaler)
      except BadDatasetForScalingException as e:
        logger.info(e)
        logger.info('Removing experiment ' + str(i) +'\n' + '='*80 + '\n')
        del experiments[i - offset]
        del reflections[i - offset]
        offset += 1
    assert len(experiments) == len(single_scalers), (len(experiments), len(single_scalers))
    assert len(experiments) == len(reflections), (len(experiments), len(reflections))
    return MultiScaler(params, experiments, single_scalers)

  @classmethod
  def create_from_targetscaler(cls, targetscaler):
    '''method to pass scalers from TargetScaler to a MultiScaler'''
    single_scalers = []
    for scaler in targetscaler.unscaled_scalers:
      #scaler.select_reflections_for_scaling(for_multi=True)
      single_scalers.append(scaler)
    single_scalers.extend(targetscaler.single_scalers)
    return MultiScaler(targetscaler.params, [targetscaler.experiments], single_scalers)

class TargetScalerFactory(object):
  'Factory for creating a targeted scaler for multiple datasets'
  @classmethod
  def create(cls, params, experiments, reflections, is_scaled_list):
    '''sort scaled and unscaled datasets to pass to TargetScaler'''
    scaled_experiments = []
    scaled_scalers = []
    unscaled_scalers = []
    offset = 0
    for i in range(len(reflections)):
      if is_scaled_list[i] is True:
        if params.scaling_options.target_model or params.scaling_options.target_mtz:
          scaled_experiments.append(experiments[i-offset])
          scaled_scalers.append(NullScalerFactory.create(params, experiments[i-offset],
            reflections[i-offset]))
        else:
          try:
            scaled_scalers.append(SingleScalerFactory.create(params,
              experiments[i-offset], reflections[i-offset], for_multi=True))
            scaled_experiments.append(experiments[i-offset])
          except BadDatasetForScalingException as e:
            logger.info(e)
            logger.info('Removing experiment ' + str(i) +'\n' + '='*80 + '\n')
            del experiments[i - offset]
            del reflections[i - offset]
            offset += 1
      else:
        try:
          unscaled_scalers.append(SingleScalerFactory.create(params,
            experiments[i-offset], reflections[i-offset], for_multi=True))
        except BadDatasetForScalingException as e:
          logger.info(e)
          logger.info('Removing experiment ' + str(i) +'\n' + '='*80 + '\n')
          del experiments[i - offset]
          del reflections[i - offset]
          offset += 1
    assert len(experiments) == len(scaled_scalers) + len(unscaled_scalers), (
      len(experiments), str(len(scaled_scalers)) + ' + ' + str(len(unscaled_scalers)))
    assert len(experiments) == len(reflections), (len(experiments), len(reflections))
    return TargetScaler(params, scaled_experiments, scaled_scalers,
      unscaled_scalers)
