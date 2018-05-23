'''
Collection of factories for creating the scalers.
'''
import logging
from libtbx.utils import Sorry
from dials.array_family import flex
from dials.algorithms.scaling.scaler import MultiScaler, TargetScaler,\
  SingleScalerBase
from dials.algorithms.scaling.scaling_utilities import calc_normE2
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

class SingleScalerFactory(object):
  """Factory for creating a scaler for a single dataset."""

  @classmethod
  def create(cls, params, experiment, reflection_table, scaled_id=0, for_multi=False):
    """Perform reflection_table preprocessing and create a SingleScaler."""

    reflection_table['id'] = flex.int(reflection_table.size(), scaled_id)
    if params.scaling_options.verbosity > 1:
      logger.info(('Preprocessing data for scaling. The id assigned to this \n'
        'dataset is {0}, and the scaling model type being applied is {1}. \n'
        ).format(scaled_id, experiment.scaling_model.id_))

    reflection_table = cls.filter_bad_reflections(reflection_table)
    if params.scaling_options.verbosity > 1:
      logger.info('%s reflections not suitable for scaling (low partiality,\n'
        'not integrated etc).\n', reflection_table.get_flags(
        reflection_table.flags.excluded_for_scaling).count(True))

    reflection_table = cls.select_optimal_intensities(reflection_table, params)

    if not 'inverse_scale_factor' in reflection_table:
      reflection_table['inverse_scale_factor'] = flex.double(
        reflection_table.size(), 1.0)
      print('setting inverse scale factors to zero')

    reflection_table = cls.filter_outliers(reflection_table, experiment,
      params)

    return SingleScalerBase(params, experiment, reflection_table, for_multi)

  @classmethod
  def filter_outliers(cls, reflections, experiment, params):
    """Calculate normalised E2 values and perform outlier rejection."""
    reflections = calc_normE2(reflections, experiment)
    if params.scaling_options.outlier_rejection != '0':
      reflections = reject_outliers(reflections,
        experiment.crystal.get_space_group(),
        params.scaling_options.outlier_rejection,
        params.scaling_options.outlier_zmax)
    return reflections

  @classmethod
  def filter_bad_reflections(cls, reflections):
    """Initial filter to select integrated reflections."""
    mask = ~reflections.get_flags(reflections.flags.integrated)
    d_mask = reflections['d'] <= 0.0
    partials_mask = reflections['partiality'] < 0.6
    reflections.set_flags(mask | partials_mask | d_mask,
      reflections.flags.excluded_for_scaling)
    return reflections

  @classmethod
  def select_optimal_intensities(cls, reflection_table, params):
    """Choose which intensities to use for scaling."""
    integration_method = params.scaling_options.integration_method
    if integration_method == 'sum' or integration_method == 'prf':
      intstr = integration_method
      conversion = flex.double(reflection_table.size(), 1.0)
      if 'partiality' in reflection_table:
        inverse_partiality = flex.double(reflection_table.size(), 1.0)
        nonzero_partiality_sel = reflection_table['partiality'] > 0.0
        good_refl = reflection_table.select(reflection_table['partiality'] > 0.0)
        inverse_partiality.set_selected(nonzero_partiality_sel.iselection(),
          1.0/good_refl['partiality'])
        conversion *= inverse_partiality
      if 'lp' in reflection_table:
        conversion *= reflection_table['lp']
      if 'qe' in reflection_table:
        conversion /= reflection_table['qe']
      elif 'dqe' in reflection_table:
        conversion /= reflection_table['dqe']
      reflection_table['intensity'] = (
        reflection_table['intensity.'+intstr+'.value'] * conversion)
      reflection_table['variance'] = (
        reflection_table['intensity.'+intstr+'.variance'] * conversion
          *conversion)
      if params.scaling_options.verbosity > 1:
        logger.info(('{0} intensity values will be used for scaling (and mtz \n'
        'output if applicable). \n').format('Profile fitted' if intstr == 'prf'
        else 'Summation integrated'))
    #perform a combined prf/sum in a similar fashion to aimless
    else:
      if params.scaling_options.verbosity > 1:
        logger.info('Intensity selection choice does not have an implementation \n'
          'using prf values.')
      params.scaling_options.integration_method = 'prf'
      return cls.select_optimal_intensities(reflection_table, params)
      '''int_prf = reflection_table['intensity.prf.value'] * conversion
      int_sum = reflection_table['intensity.sum.value'] * conversion
      var_prf = reflection_table['intensity.prf.variance'] * (conversion**2)
      var_sum = reflection_table['intensity.sum.variance'] * (conversion**2)
      Imid = max(int_sum)/2.0
      weight = 1.0/(1.0 + ((int_prf/Imid)**3))
      reflection_table['intensity'] = ((weight * int_prf)
        + ((1.0 - weight) * int_sum))
      reflection_table['variance'] = ((weight * var_prf)
        + ((1.0 - weight) * var_sum))
      msg = ('Combined profile/summation intensity values will be used for {sep}'
      'scaling, with an Imid of {0}. {sep}').format(Imid, sep='\n')
      logger.info(msg)'''
    variance_mask = reflection_table['variance'] <= 0.0
    reflection_table.set_flags(variance_mask,
      reflection_table.flags.excluded_for_scaling)
    return reflection_table

class NullScalerFactory(object):
  'Factory for creating null scaler'
  @classmethod
  def create(cls, params, experiment, reflection, scaled_id=0):
    """Return Null Scaler."""
    from dials.algorithms.scaling.scaler import NullScaler
    return NullScaler(params, experiment, reflection, scaled_id)

class MultiScalerFactory(object):
  'Factory for creating a scaler for multiple datasets'
  @classmethod
  def create(cls, params, experiments, reflections):
    '''create a list of single scalers to pass to a MultiScaler. For scaled_id,
    we just need unique values, not necessarily the same as previously.'''
    single_scalers = []
    for i, (reflection, experiment) in enumerate(zip(reflections, experiments)):
      single_scalers.append(SingleScalerFactory.create(
        params, experiment, reflection, scaled_id=i, for_multi=True))
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
  def create(cls, params, experiments, reflections, is_scaled_list):
    '''sort scaled and unscaled datasets to pass to TargetScaler'''
    scaled_experiments = []
    scaled_scalers = []
    unscaled_scalers = []
    for i, (experiment, reflection) in enumerate(zip(experiments, reflections)):
      if is_scaled_list[i] is True:
        if params.scaling_options.target_model:
          scaled_experiments.append(experiment)
          scaled_scalers.append(NullScalerFactory.create(params, experiment,
            reflection, scaled_id=i))
        else:
          scaled_experiments.append(experiment)
          scaled_scalers.append(SingleScalerFactory.create(params, experiment,
            reflection, scaled_id=i, for_multi=True))
      else:
        unscaled_scalers.append(SingleScalerFactory.create(params, experiment,
          reflection, scaled_id=i, for_multi=True))
    return TargetScaler(params, scaled_experiments, scaled_scalers,
      unscaled_scalers)
