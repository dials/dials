'''
Collection of factories for creating the scalers.
'''
import logging
from libtbx.utils import Sorry
from dials.array_family import flex
from dials.algorithms.scaling.scaler import MultiScaler, TargetScaler,\
  SingleScalerBase
from dials.algorithms.scaling.scaling_utilities import quasi_normalisation
from dials.algorithms.scaling.scaling_library import choose_scaling_intensities
from dials.algorithms.scaling.outlier_rejection import reject_outliers
logger = logging.getLogger('dials')

def create_scaler(params, experiments, reflections):
  """Read an experimentlist and list of reflection tables and return
    an appropriate scaler."""
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
    mask = ~reflections.get_flags(reflections.flags.integrated, all=False)
    if 'd' in reflections:
      d_mask = reflections['d'] <= 0.0
      mask = mask | d_mask
    reflections.set_flags(mask, reflections.flags.excluded_for_scaling)
    return reflections

class SingleScalerFactory(ScalerFactory):
  """Factory for creating a scaler for a single dataset."""

  @classmethod
  def create(cls, params, experiment, reflection_table, scaled_id=0, for_multi=False):
    """Perform reflection_table preprocessing and create a SingleScaler."""

    if experiment.identifier:
      assert experiment.identifier == \
        reflection_table.experiment_identifiers().values()[0]
      if params.scaling_options.verbosity > 1:
        logger.info('The experiment identifier for this dataset is %s',
          experiment.identifier)
    else:
      reflection_table['id'] = flex.int(reflection_table.size(), scaled_id)
    if params.scaling_options.verbosity > 1:
      logger.info('Preprocessing data for scaling. The id assigned to this \n'
        'dataset is %s, and the scaling model type being applied is %s. \n' %
        (reflection_table['id'][0], experiment.scaling_model.id_))

    reflection_table = cls.filter_bad_reflections(reflection_table)
    if params.scaling_options.verbosity > 1:
      logger.info('%s reflections not suitable for scaling (bad d-value,\n'
        'not integrated etc).\n', reflection_table.get_flags(
        reflection_table.flags.excluded_for_scaling).count(True))

    if not 'inverse_scale_factor' in reflection_table:
      reflection_table['inverse_scale_factor'] = flex.double(
        reflection_table.size(), 1.0)

    reflection_table = choose_scaling_intensities(reflection_table,
      params.reflection_selection.intensity_choice)

    reflection_table = cls.filter_outliers(reflection_table, experiment,
      params)

    return SingleScalerBase(params, experiment, reflection_table, for_multi)

  @classmethod
  def filter_outliers(cls, reflections, experiment, params):
    """Calculate normalised E2 values and perform outlier rejection."""
    reflections = quasi_normalisation(reflections, experiment)
    if params.scaling_options.outlier_rejection:
      reflections = reject_outliers([reflections],
        experiment.crystal.get_space_group(),
        params.scaling_options.outlier_rejection,
        params.scaling_options.outlier_zmax)[0]
    return reflections

class NullScalerFactory(ScalerFactory):
  'Factory for creating null scaler'
  @classmethod
  def create(cls, params, experiment, reflection, scaled_id=0):
    """Return Null Scaler."""
    from dials.algorithms.scaling.scaler import NullScaler
    logger.info('Preprocessing target dataset for scaling. \n')
    reflection = cls.filter_bad_reflections(reflection)
    variance_mask = reflection['variance'] <= 0.0
    reflection.set_flags(variance_mask, reflection.flags.excluded_for_scaling)
    logger.info('%s reflections not suitable for scaling (bad d-value,\n'
        'not integrated etc).\n', reflection.get_flags(
        reflection.flags.excluded_for_scaling).count(True))
    return NullScaler(params, experiment, reflection, scaled_id)

class MultiScalerFactory(object):
  'Factory for creating a scaler for multiple datasets'
  @classmethod
  def create(cls, params, experiments, reflections, dataset_ids=None):
    '''create a list of single scalers to pass to a MultiScaler. For scaled_id,
    we just need unique values, not necessarily the same as previously.'''
    single_scalers = []
    if not dataset_ids:
      dataset_ids = range(len(reflections))
    for (data_id, reflection, experiment) in zip(dataset_ids, reflections, experiments):
      single_scalers.append(SingleScalerFactory.create(
        params, experiment, reflection, scaled_id=data_id, for_multi=True))
    return MultiScaler(params, experiments, single_scalers)

  @classmethod
  def create_from_targetscaler(cls, targetscaler):
    '''method to pass scalers from TargetScaler to a MultiScaler'''
    single_scalers = targetscaler.single_scalers
    for scaler in targetscaler.unscaled_scalers:
      scaler.select_reflections_for_scaling(for_multi=True)
      single_scalers.append(scaler)
    return MultiScaler(targetscaler.params, [targetscaler.experiments], single_scalers)

class TargetScalerFactory(object):
  'Factory for creating a targeted scaler for multiple datasets'
  @classmethod
  def create(cls, params, experiments, reflections, is_scaled_list, dataset_ids=None):
    '''sort scaled and unscaled datasets to pass to TargetScaler'''
    scaled_experiments = []
    scaled_scalers = []
    unscaled_scalers = []
    if not dataset_ids:
      dataset_ids = range(len(reflections))
    for i, (experiment, reflection) in enumerate(zip(experiments, reflections)):
      if is_scaled_list[i] is True:
        if params.scaling_options.target_model or params.scaling_options.target_mtz:
          scaled_experiments.append(experiment)
          scaled_scalers.append(NullScalerFactory.create(params, experiment,
            reflection, scaled_id=dataset_ids[i]))
        else:
          scaled_experiments.append(experiment)
          scaled_scalers.append(SingleScalerFactory.create(params, experiment,
            reflection, scaled_id=dataset_ids[i], for_multi=True))
      else:
        unscaled_scalers.append(SingleScalerFactory.create(params, experiment,
          reflection, scaled_id=dataset_ids[i], for_multi=True))
    return TargetScaler(params, scaled_experiments, scaled_scalers,
      unscaled_scalers)
