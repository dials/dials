'''
Collection of factories for creating the scalers.
'''
import logging
import numpy as np
from dials.array_family import flex
from cctbx import miller, crystal
from scitbx import sparse
import iotbx.merging_statistics
from dials.algorithms.scaling.target_function import \
  target_function, target_function_fixedIh
from dials.algorithms.scaling.basis_functions import basis_function
from dials.algorithms.scaling.scaling_utilities import calc_s2d, sph_harm_table
from dials.algorithms.scaling.Wilson_outlier_test import (
  calculate_wilson_outliers, calc_normE2)
from dials.algorithms.scaling.reflection_weighting import Weighting
from dials.algorithms.scaling.target_Ih import SingleIhTable, JointIhTable, IhTableBase
from dials.algorithms.scaling.minimiser_functions import error_scale_LBFGSoptimiser
from dials.algorithms.scaling.aimless_outlier_rejection import reject_outliers
from dxtbx.model import Crystal

logger = logging.getLogger('dials')

class Factory(object):
  '''
  Factory for creating Scalers.
  '''
  @classmethod
  def create(cls, params, experiments, reflections):
    '''
    create the scaling model defined by the params.
    '''
    if len(reflections) == 1:
      scaler = SingleScalerFactory.create(params, experiments[0], reflections[0])
    else:
      is_scaled_list = cls.is_scaled(reflections)
      #check if we want to do targeted scaling or normal scaling
      n_scaled = is_scaled_list.count(True)
      if (params.scaling_options.target is True and n_scaled > 0
        and n_scaled < len(reflections)):
        scaler = TargetScalerFactory.create(params, experiments, reflections,
          is_scaled_list)
      elif len(reflections) > 1:
        scaler = MultiScalerFactory.create(params, experiments, reflections)
      else:
        assert 0, 'no reflection tables found to create the scaler'
    return scaler

  @classmethod
  def is_scaled(cls, reflections):
    '''inspect reflection table to see if it already has scale factors.'''
    is_already_scaled = []
    for reflection_table in reflections:
      if 'scaled_id' in reflection_table.keys():
        is_already_scaled.append(True)
      else:
        is_already_scaled.append(False)
    return is_already_scaled


class SingleScalerFactory(object):
  'Factory for creating a scaler for a single dataset'
  @classmethod
  def create(cls, params, experiment, reflection, scaled_id=0):
    '''create a single scaler with the relevant parameterisation'''
    if experiment.scaling_model.id_ == "aimless":
      return AimlessScaler(params, experiment, reflection, scaled_id)
    elif experiment.scaling_model.id_ == "KB":
      return KBScaler(params, experiment, reflection, scaled_id)
    elif experiment.scaling_model.id_ == "xscale":
      return XscaleScaler(params, experiment, reflection, scaled_id)
    else:
      assert 0, "Scaling model __id__ not recognised."

class MultiScalerFactory(object):
  'Factory for creating a scaler for multiple datasets'
  @classmethod
  def create(cls, params, experiments, reflections):
    '''create a list of single scalers to pass to a MultiScaler'''
    single_scalers = []
    for i, (reflection, experiment) in enumerate(zip(reflections, experiments)):
      single_scalers.append(SingleScalerFactory.create(
        params, experiment, reflection, scaled_id=i))
    return MultiScaler(params, experiments, single_scalers)

  @classmethod
  def create_from_targetscaler(cls, targetscaler):
    '''pass a TargetScaler to a MultiScaler'''
    single_scalers = targetscaler.single_scalers
    for scaler in targetscaler.unscaled_scalers:
      single_scalers.append(scaler)
    return MultiScaler(targetscaler.params, [targetscaler.experiments], single_scalers)

class TargetScalerFactory(object):
  'Factory for creating a targeted scaler for multiple datasets'
  @classmethod
  def create(cls, params, experiments, reflections, is_scaled_list):
    '''sort scaled and unscaled datasets to pass to TargetScaler'''
    scaled_experiments = []
    scaled_scalers = []
    unscaled_experiments = []
    unscaled_scalers = []
    for i, reflection in enumerate(reflections):
      if is_scaled_list[i] is True:
        scaled_experiments.append(experiments[i])
        scaled_scalers.append(SingleScalerFactory.create(params, experiments[i],
          reflection, scaled_id=i))
      else:
        unscaled_experiments.append(experiments[i])
        unscaled_scalers.append(SingleScalerFactory.create(params, experiments[i],
          reflection, scaled_id=i))
    #if len(unscaled_scalers) > 1:
    #  assert 0, """method for targeted scaling of multiple datasets
    #  against a target not yet implemented"""
    return TargetScaler(params, scaled_experiments, scaled_scalers,
      unscaled_experiments, unscaled_scalers)


class ScalerUtilities(object):
  '''Base class for all Scalers (single and multiple)'''
  def __init__(self):
    'General attributes relevant for all parameterisations'
    self._experiments = None
    self._params = None
    self._reflection_table = []
    self._outlier_table = flex.reflection_table()
    self._Ih_table = None
    self._initial_keys = []
    self._g_parameterisation = {}

  @property
  def Ih_table(self):
    return self._Ih_table

  @Ih_table.setter
  def Ih_table(self, new_Ih_table):
    assert isinstance(new_Ih_table, IhTableBase)
    self._Ih_table = new_Ih_table

  @property
  def g_parameterisation(self):
    '''dict to hold {param_name:SF} object pairs to pass to apm'''
    return self._g_parameterisation

  @property
  def experiments(self):
    return self._experiments

  @property
  def reflection_table(self):
    return self._reflection_table

  @property
  def outlier_table(self):
    return self._outlier_table

  @property
  def params(self):
    return self._params

  @property
  def initial_keys(self):
    '''list of initial reflection table keys.'''
    return self._initial_keys

  def clean_reflection_table(self):
    '''remove additional added columns that are not required for output'''
    self._initial_keys.append('inverse_scale_factor')
    self._initial_keys.append('Ih_values')
    for key in self.reflection_table.keys():
      if not key in self._initial_keys:
        del self._reflection_table[key]

  def update_for_minimisation(self, apm):
    '''update the scale factors and Ih for the next iteration of minimisation'''
    basis_fn = self.get_basis_function(apm)
    apm.active_derivatives = basis_fn[1]
    self.Ih_table.inverse_scale_factors = basis_fn[0]
    self.Ih_table.calc_Ih()

  def get_target_function(self, apm):
    '''call the target function'''
    return target_function(self, apm).return_targets()

  def get_basis_function(self, apm):
    '''call thebasis function'''
    return basis_function(self, apm).return_basis()

  def save_reflection_table(self, filename):
    ''' Save the reflections to file. '''
    self.reflection_table.as_pickle(filename)

  def save_outlier_table(self, filename):
    ''' Save the reflections to file. '''
    self.outlier_table.as_pickle(filename)

  @classmethod
  def _scaling_subset(cls, reflection_table, params, error_model_params=None):
    '''select the reflections with non-zero weight and update scale weights
    object.'''
    weights_for_scaling = cls._update_weights_for_scaling(reflection_table,
      params, error_model_params=error_model_params)
    sel = weights_for_scaling.weights > 0.0
    sel1 = reflection_table['Esq'] > params.reflection_selection.E2min
    sel2 = reflection_table['Esq'] < params.reflection_selection.E2max
    selection = sel & sel1 & sel2
    reflections_for_scaling = reflection_table.select(selection)
    weights_for_scaling.weights = weights_for_scaling.weights.select(selection)
    msg = ('{0} reflections were selected for scale factor determination {sep}'
      'out of {5} reflections. This was based on selection criteria of {sep}'
      'E2min = {1}, E2max = {2}, Isigma_min = {3}, dmin = {4}. {sep}').format(
      reflections_for_scaling.size(), params.reflection_selection.E2min,
      params.reflection_selection.E2max, params.reflection_selection.Isigma_min,
      params.reflection_selection.d_min, len(reflection_table), sep='\n')
    logger.info(msg)
    return reflections_for_scaling, weights_for_scaling, selection

  @staticmethod
  def _update_weights_for_scaling(reflection_table, params,
    weights_filter=True, error_model_params=None):
    '''set the weights of each reflection to be used in scaling'''
    weights_for_scaling = Weighting(reflection_table)
    logger.info('Updating the weights associated with the intensities. \n')
    if weights_filter:
      weights_for_scaling.apply_Isigma_cutoff(reflection_table,
        params.reflection_selection.Isigma_min)
      weights_for_scaling.apply_dmin_cutoff(reflection_table,
        params.reflection_selection.d_min)
    weights_for_scaling.remove_wilson_outliers(reflection_table)
    #if params.weighting.tukey_biweighting and Ih_table:
    #  weights_for_scaling.tukey_biweighting(Ih_table)
    if error_model_params:
      weights_for_scaling.apply_aimless_error_model(reflection_table,
        error_model_params)
    return weights_for_scaling

  def calc_merging_statistics(self):
    u_c = self.experiments.crystal.get_unit_cell().parameters()
    if self.params.scaling_options.force_space_group:
      s_g_symbol = self.params.scaling_options.force_space_group
      crystal_symmetry = crystal.symmetry(unit_cell=u_c,
        space_group_symbol=s_g_symbol)
    else:
      s_g = self.experiments.crystal.get_space_group()
      crystal_symmetry = crystal.symmetry(unit_cell=u_c,
        space_group=s_g)
    miller_set = miller.set(crystal_symmetry=crystal_symmetry,
      indices=self.Ih_table.miller_index, anomalous_flag=False)
    scaled_intensities = self.Ih_table.intensities/self.Ih_table.inverse_scale_factors
    sigmas = (self.Ih_table.variances**0.5)/self.Ih_table.inverse_scale_factors
    i_obs = miller.array(miller_set, data=scaled_intensities, sigmas=sigmas)
    i_obs.set_observation_type_xray_intensity()
    scaled_ids = list(set(self.reflection_table['scaled_id']))
    result = iotbx.merging_statistics.dataset_statistics(
      i_obs=i_obs, n_bins=20, anomalous=False, sigma_filtering=None,
      use_internal_variance=True)
      #eliminate_sys_absent=False, sigma_filtering=None)
    return ([result], scaled_ids[0])

  def calc_correlation(self):
    if isinstance(self, MultiScaler):
      array_list = []
      for scaler in self.single_scalers:
        u_c = scaler.experiments.crystal.get_unit_cell().parameters()
        if self.params.scaling_options.force_space_group:
          s_g_symbol = self.params.scaling_options.force_space_group
          crystal_symmetry = crystal.symmetry(unit_cell=u_c,
            space_group_symbol=s_g_symbol)
        else:
          s_g = scaler.experiments.crystal.get_space_group()
          crystal_symmetry = crystal.symmetry(unit_cell=u_c,
            space_group=s_g)
        miller_set = miller.set(crystal_symmetry=crystal_symmetry,
          indices=scaler.Ih_table.asu_miller_index, anomalous_flag=False)
        scaled_intensities = scaler.Ih_table.intensities/scaler.Ih_table.inverse_scale_factors
        sigmas = (scaler.Ih_table.variances**0.5)/scaler.Ih_table.inverse_scale_factors
        i_obs = miller.array(miller_set, data=scaled_intensities, sigmas=sigmas)
        array_list.append(i_obs)
      correl_list = []
      for array_1 in array_list:
        for array_2 in array_list:
          correl = array_1.correlation(other=array_2)
          correl_list.append(correl.coefficient())
      return correl_list
    return None



class SingleScaler(ScalerUtilities):
  '''
  Parent class for single dataset Scalers, containing a standard
  setup routine for the reflection_table - takes in params, experiment
  and reflection.
  '''
  def __init__(self, params, experiment, reflection, scaled_id=0):
    logger.info('\nInitialising a Single Scaler instance. \n')
    super(SingleScaler, self).__init__()
    self._experiments = experiment
    self._corrections = experiment.scaling_model.configdict['corrections']
    self._params = params
    logger.info("Dataset id for this reflection table is %s." % scaled_id)
    logger.info(('The type of scaling model being applied to this dataset {sep}'
      'is {0}. {sep}').format(self.experiments.scaling_model.id_, sep='\n'))
    reflection['scaled_id'] = flex.int([scaled_id]*len(reflection))
    self._initial_keys = [key for key in reflection.keys()]
    #choose intensities, map to asu, assign unique refl. index
    reflection_table = self._reflection_table_setup(self._initial_keys, reflection)
    reflection_table = self._select_optimal_intensities(reflection_table, self.params)
    reflection_table = self._map_indices_to_asu(reflection_table,
      self.experiments, self.params)
    #calculate values for later filtering, but don't filter here!!!
    reflection_table = calc_normE2(reflection_table, self.experiments)
    reflection_table['wilson_outlier_flag'] = calculate_wilson_outliers(
      reflection_table)
    self._reflection_table = reflection_table
    

  @property
  def corrections(self):
    return self._corrections

  @staticmethod
  def _reflection_table_setup(initial_keys, reflections):
    'initial filter to select integrated reflections'
    reflections = reflections.select(reflections.get_flags(
      reflections.flags.integrated))
    if 'intensity.prf.variance' in initial_keys:
      reflections = reflections.select(reflections['intensity.prf.variance'] > 0)
    if 'intensity.sum.variance' in initial_keys:
      reflections = reflections.select(reflections['intensity.sum.variance'] > 0)
    if not 'inverse_scale_factor' in initial_keys:
      reflections['inverse_scale_factor'] = (flex.double([1.0] * len(reflections)))
    reflections = reflections.select(reflections['partiality'] > 0.95)
    return reflections

  @staticmethod
  def _map_indices_to_asu(reflection_table, experiments, params):
    '''Create a miller_set object, map to the asu and create a sorted
       reflection table, sorted by asu miller index'''
    u_c = experiments.crystal.get_unit_cell().parameters()
    if params.scaling_options.force_space_group:
      sg_from_file = experiments.crystal.get_space_group().info()
      s_g_symbol = params.scaling_options.force_space_group
      crystal_symmetry = crystal.symmetry(unit_cell=u_c,
        space_group_symbol=s_g_symbol)
      msg = ('WARNING: Manually overriding space group from {0} to {1}. {sep}'
        'If the reflection indexing in these space groups is different, {sep}'
        'bad things may happen!!! {sep}').format(sg_from_file, s_g_symbol, sep='\n')
      logger.info(msg)
    else:
      s_g = experiments.crystal.get_space_group()
      crystal_symmetry = crystal.symmetry(unit_cell=u_c, space_group=s_g)
    miller_set = miller.set(crystal_symmetry=crystal_symmetry,
      indices=reflection_table['miller_index'], anomalous_flag=False)
    miller_set_in_asu = miller_set.map_to_asu()
    reflection_table["asu_miller_index"] = miller_set_in_asu.indices()
    permuted = (miller_set.map_to_asu()).sort_permutation(by_value='packed_indices')
    reflection_table = reflection_table.select(permuted)
    return reflection_table

  @classmethod
  def _select_optimal_intensities(cls, reflection_table, params):
    '''method to choose which intensities to use for scaling'''
    if (params.scaling_options.integration_method == 'sum' or
        params.scaling_options.integration_method == 'prf'):
      intstr = params.scaling_options.integration_method
      reflection_table['intensity'] = (reflection_table['intensity.'+intstr+'.value']
                                       * reflection_table['lp']
                                       / reflection_table['dqe'])
      reflection_table['variance'] = (reflection_table['intensity.'+intstr+'.variance']
                                      * (reflection_table['lp']**2)
                                      / (reflection_table['dqe']**2))
      logger.info(('{0} intensity values will be used for scaling. {sep}').format(
        'Profile fitted' if intstr == 'prf' else 'Summation integrated', sep='\n'))
    #perform a combined prf/sum in a similar fashion to aimless
    elif params.scaling_options.integration_method == 'combine':
      int_prf = (reflection_table['intensity.prf.value']
                 * reflection_table['lp'] / reflection_table['dqe'])
      int_sum = (reflection_table['intensity.sum.value']
                 * reflection_table['lp'] / reflection_table['dqe'])
      var_prf = (reflection_table['intensity.prf.variance']
                 * (reflection_table['lp']**2) / (reflection_table['dqe']**2))
      var_sum = (reflection_table['intensity.sum.variance']
                 * (reflection_table['lp']**2) / (reflection_table['dqe']**2))
      Imid = max(int_sum)/2.0
      weight = 1.0/(1.0 + ((int_prf/Imid)**3))
      reflection_table['intensity'] = ((weight * int_prf) + ((1.0 - weight) * int_sum))
      reflection_table['variance'] = ((weight * var_prf) + ((1.0 - weight) * var_sum))
      msg = ('Combined profile/summation intensity values will be used for {sep}'
      'scaling, with an Imid of {0}. {sep}').format(Imid, sep='\n')
      logger.info(msg)
    else:
      logger.info('Invalid integration_method choice, using default profile fitted intensities')
      params.scaling_options.integration_method = 'prf'
      cls._select_optimal_intensities(reflection_table, params)
    return reflection_table

  def expand_scales_to_all_reflections(self, caller=None):
    self._reflection_table['inverse_scale_factor'] = self.calc_expanded_scales()
    logger.info(('Scale factors determined during minimisation have now been applied {sep}'
      'to all reflections. {sep}').format(sep='\n'))
    weights = self._update_weights_for_scaling(self.reflection_table,
      self.params, weights_filter=False, error_model_params=None).weights
    #now create an Ih table that doesn't include any outliers
    self._Ih_table = SingleIhTable(self.reflection_table, weights)
    outlier_sel = weights == 0.0
    self._outlier_table.extend(self.reflection_table.select(outlier_sel))
    self.Ih_table = self.Ih_table.select(self.Ih_table.weights != 0.0)
    self._reflection_table = self._reflection_table.select(weights != 0.0)
    #if self.params.weighting.tukey_biweighting:
    #  self.Ih_table.apply_tukey_biweighting()
    #  self.Ih_table.calc_Ih()
    self._reflection_table['Ih_values'] = self.Ih_table.Ih_values
    if (self.params.scaling_options.reject_outliers and
      not isinstance(caller, MultiScaler)):
      self.round_of_outlier_rejection()
    logger.info('A new best estimate for I_h for all reflections has now been calculated. \n')

  def update_error_model(self):
    '''apply a correction to try to improve the error estimate.'''
    self.params.weighting.error_model_params = (
      error_scale_LBFGSoptimiser(self.Ih_table, flex.double([1.0, 0.01])).x)
    self.Ih_table.update_aimless_error_model(self.params.weighting.error_model_params)

  def round_of_outlier_rejection(self):
    '''calculate outliers from the reflections in the Ih_table,
    and use these to filter the reflection table and Ih_table.'''
    sel = reject_outliers(self, self.params.scaling_options.outlier_zmax)
    outliers_sel = ~sel
    self._outlier_table.extend(self._reflection_table.select(outliers_sel))
    self._reflection_table = self._reflection_table.select(sel)
    self.Ih_table = self.Ih_table.select(sel)
    self._reflection_table['Ih_values'] = self.Ih_table.Ih_values
    msg = ('A round of outlier rejection has been performed, {0} outliers {sep}'
        'were found which have been removed from the dataset. {sep}'.format(
        sel.count(False), sep='\n'))
    logger.info(msg)

  def apply_selection_to_SFs(self, sel):
    '''Updates data from within current SF objects using the given selection.
       Required for targeted scaling. To be filled in by subclasses'''
    pass

  def calc_expanded_scales(self):
    '''calculate the scale factors for all reflections from the model.
    To be filled in by subclasses.'''
    pass

  def normalise_scale_component(self):
    pass

  def normalise_decay_component(self):
    pass

  def calc_absorption_constraint(self, apm):
    '''calculates a constraint for the spherical harmonic absorption correction.
    Should only be called from target function if g_absorption in active params.'''
    if self.id_ == 'aimless':
      if 'g_absorption' in apm.active_parameterisation:
        idx = apm.active_parameterisation.index('g_absorption')
        start_idx = apm.n_cumul_params_list[idx]
        end_idx = apm.n_cumul_params_list[idx+1]
        weight = self.params.parameterisation.surface_weight
        abs_params = apm.x[start_idx:end_idx]
        residual = (weight * (abs_params)**2)
        gradient = (2.0 * weight * abs_params)
        #return a gradient vector to be added to that calculated in target function
        gradient_vector = flex.double([])
        for i, param in enumerate(apm.active_parameterisation):
          if param != 'g_absorption':
            gradient_vector.extend(flex.double([0.0]*apm.n_params_list[i]))
          elif param == 'g_absorption':
            gradient_vector.extend(gradient)
        return (residual, gradient_vector)
      else:
        return (flex.double([0.0]), flex.double([0.0]*apm.n_active_params))
    return (flex.double([0.0]), flex.double([0.0]*apm.n_active_params))

class KBScaler(SingleScaler):
  '''
  Scaler for single dataset using simple KB parameterisation.
  '''

  id_ = 'KB'

  def __init__(self, params, experiment, reflection, scaled_id=0):
    super(KBScaler, self).__init__(params, experiment, reflection, scaled_id)
    (self.g_scale, self.g_decay) = (None, None)
    self._initialise_scale_factors()
    self._select_reflections_for_scaling()
    logger.info('Completed initialisation of KB Scaler. \n' + '*'*40 + '\n')

  def _initialise_scale_factors(self):
    if 'scale' in self.corrections:
      self.g_scale = self.experiments.scaling_model.components['scale']
      self._g_parameterisation['g_scale'] = self.g_scale
    if 'decay' in self.corrections:
      self.g_decay = self.experiments.scaling_model.components['decay']
      self._g_parameterisation['g_decay'] = self.g_decay

  def _select_reflections_for_scaling(self):
    (reflections_for_scaling, weights_for_scaling, _) = (
      self._scaling_subset(self.reflection_table, self.params))
    self._Ih_table = SingleIhTable(reflections_for_scaling, weights_for_scaling.weights)
    if 'scale' in self.corrections:
      self.g_scale.update_reflection_data(n_refl=self.Ih_table.size)
    if 'decay' in self.corrections:
      self.g_decay.update_reflection_data(dvalues=reflections_for_scaling['d'])

  def apply_selection_to_SFs(self, sel):
    '''Updates data from within current SF objects using the given selection.
       Required for targeted scaling.'''
    if 'decay' in self.corrections:
      self.g_decay.update_reflection_data(dvalues=self.g_decay.d_values.select(sel))
    if 'scale' in self.corrections:
      self.g_scale.update_reflection_data(n_refl=sel.count(True))

  def calc_expanded_scales(self):
    '''calculate'''
    expanded_scale_factors = flex.double([1.0]*len(self.reflection_table))
    if 'scale' in self.corrections:
      self.g_scale.update_reflection_data(n_refl=len(self.reflection_table))
      self.g_scale.calculate_scales_and_derivatives()
      expanded_scale_factors *= self.g_scale.inverse_scales
      logger.info(('Scale factor K = {0:.4f} determined during minimisation {sep}'
        'has now been applied to all reflections. {sep}').format(
        list(self.g_scale.parameters)[0], sep='\n'))
    if 'decay' in self.corrections:
      self.g_decay.update_reflection_data(dvalues=self.reflection_table['d'])
      self.g_decay.calculate_scales_and_derivatives()
      expanded_scale_factors *= self.g_decay.inverse_scales
      logger.info(('B-factor B = {0:.4f} determined during minimisation {sep}'
        'has now been applied to all reflections. {sep}').format(
        list(self.g_decay.parameters)[0], sep='\n'))
    return expanded_scale_factors


class AimlessScaler(SingleScaler):
  '''
  Scaler for single dataset using aimless-like parameterisation.
  '''

  id_ = 'aimless'

  def __init__(self, params, experiment, reflection, scaled_id=0):
    super(AimlessScaler, self).__init__(params, experiment, reflection, scaled_id)
    (self.g_absorption, self.g_scale, self.g_decay) = (None, None, None)
    self.sph_harm_table = None
    #determine outliers, initialise scalefactors and extract data for scaling
    if self.params.scaling_options.reject_outliers:
      self._Ih_table = SingleIhTable(self.reflection_table)
      self.round_of_outlier_rejection()
    self._initialise_scale_factors()
    self._select_reflections_for_scaling()
    logger.info('Completed initialisation of AimlessScaler. \n' + '*'*40 + '\n')

  def _initialise_scale_factors(self):
    '''initialise scale factors and add to self.active_parameters'''
    if 'scale' in self.corrections:
      self._initialise_scale_term(self.reflection_table)
    if 'decay' in self.corrections:
      self._initialise_decay_term(self.reflection_table)
    if 'absorption' in self.corrections:
      self._initialise_absorption_term(self.reflection_table,
        self.experiments.scaling_model.configdict['lmax'])
        #self.params.parameterisation.lmax)
    else:
      #use calc_s2d to calculate phi values.
      self._reflection_table = calc_s2d(self.reflection_table, self.experiments)

  def _initialise_scale_term(self, refl_table):
    '''calculate the 'normalised rotation angle', and initialise a SmoothScaleFactor'''
    refl_table['norm_rot_angle'] = (refl_table['xyzobs.px.value'].parts()[2]
      * self.experiments.scaling_model.scale_normalisation_factor)
    self.g_scale = self.experiments.scaling_model.components['scale']
    self._g_parameterisation['g_scale'] = self.g_scale
    msg = ('The scale term ScaleFactor object was successfully initialised. {sep}'
      '{0} parameters will be used to parameterise the time-dependent scale. {sep}'
        ).format(self.g_scale.n_params, sep='\n')
    logger.info(msg)

  def _initialise_decay_term(self, refl_table):
    '''calculate the 'normalised time', and initialise a SmoothBScaleFactor'''
    refl_table['norm_time_values'] = (refl_table['xyzobs.px.value'].parts()[2]
      * self.experiments.scaling_model.decay_normalisation_factor)
    self.g_decay = self.experiments.scaling_model.components['decay']
    self._g_parameterisation['g_decay'] = self.g_decay
    msg = ('The decay term ScaleFactor object was successfully initialised. {sep}'
      '{0} parameters will be used to parameterise the time-dependent decay. {sep}'
      ).format(self.g_decay.n_params, sep='\n')
    logger.info(msg)

  def _initialise_absorption_term(self, reflection_table, lmax):
    reflection_table = calc_s2d(reflection_table, self.experiments)
    self.sph_harm_table = sph_harm_table(reflection_table, lmax)
    self.g_absorption = self.experiments.scaling_model.components['absorption']
    self._g_parameterisation['g_absorption'] = self.g_absorption
    msg = ('The absorption term ScaleFactor object was successfully initialised. {sep}'
      'The absorption term will be parameterised by a set of spherical {sep}'
      'harmonics up to an lmax of {0} ({1} parameters). {sep}'
      ).format(lmax, self.g_absorption.n_params, sep='\n')
    logger.info(msg)

  def _select_reflections_for_scaling(self):
    (refl_for_scaling, weights_for_scaling, selection) = (
      self._scaling_subset(self.reflection_table, self.params,
      error_model_params=self.params.weighting.error_model_params))
    self._Ih_table = SingleIhTable(refl_for_scaling, weights_for_scaling.weights)
    if self.params.weighting.tukey_biweighting:
      self.Ih_table.apply_tukey_biweighting()
    '''refactor the next two operations into extract_reflections?
    reset the normalised values within the scale_factor object to current'''
    if 'scale' in self.corrections:
      self.g_scale.update_reflection_data(refl_for_scaling['norm_rot_angle'])
    if 'decay' in self.corrections:
      self.g_decay.update_reflection_data(dvalues=refl_for_scaling['d'],
        normalised_values=refl_for_scaling['norm_time_values'])
    if 'absorption' in self.corrections:
      sph_harm_table_T = self.sph_harm_table.transpose()
      sel_sph_harm_table = sph_harm_table_T.select_columns(selection.iselection())
      self.g_absorption.update_reflection_data(sel_sph_harm_table.transpose())

  def apply_selection_to_SFs(self, sel):
    '''Updates data from within current SF objects using the given selection.
       Required for targeted scaling.'''
    if 'scale' in self.corrections:
      self.g_scale.update_reflection_data(self.g_scale.normalised_values.select(sel))
    if 'decay' in self.corrections:
      self.g_decay.update_reflection_data(dvalues=self.g_decay.d_values.select(sel),
        normalised_values=self.g_decay.normalised_values.select(sel))
    if 'absorption' in self.corrections:
      sph_harm_table_T = self.g_absorption.harmonic_values.transpose()
      sel_sph_harm_table = sph_harm_table_T.select_columns(sel.iselection())
      self.g_absorption.update_reflection_data(sel_sph_harm_table.transpose())

  

  def normalise_scale_component(self):
    '''Method to do an invariant rescale of the scale at t=0 to one.'''
    sel = (self.g_scale.normalised_values == min(self.g_scale.normalised_values))
    initial_scale = self.g_scale.inverse_scales.select(sel)[0]
    self.g_scale.parameters /= initial_scale
    self.g_scale.calculate_scales_and_derivatives()
    logger.info('Rescaled the scale component so that the initial scale is 1.\n')

  def normalise_decay_component(self):
    '''Method to do an invariant rescale of the max B to zero.'''
    maxB = max(flex.double(np.log(self.g_decay.inverse_scales))
                 * 2.0 * (self.g_decay.d_values**2))
    self.g_decay.parameters -= flex.double([maxB] * self.g_decay.n_params)
    self.g_decay.calculate_scales_and_derivatives()
    logger.info('Rescaled the decay component so that the max B is 0.\n')

  def calc_expanded_scales(self):
    expanded_scale_factors = flex.double([1.0]*len(self.reflection_table))
    if 'scale' in self.corrections:
      self.g_scale.update_reflection_data(self.reflection_table['norm_rot_angle'])
      self.g_scale.calculate_scales()
      expanded_scale_factors *= self.g_scale.inverse_scales
    if 'decay' in self.corrections:
      self.g_decay.update_reflection_data(dvalues=self.reflection_table['d'],
        normalised_values=self.reflection_table['norm_time_values'])
      self.g_decay.calculate_scales()
      expanded_scale_factors *= self.g_decay.inverse_scales
    if 'absorption' in self.corrections:
      self.g_absorption.update_reflection_data(self.sph_harm_table)
      self.g_absorption.calculate_scales_and_derivatives()
      expanded_scale_factors *= self.g_absorption.inverse_scales
    return expanded_scale_factors

  def clean_reflection_table(self):
    self._initial_keys.append('phi')
    super(AimlessScaler, self).clean_reflection_table()


class XscaleScaler(SingleScaler):
  '''
  Scaler for single dataset using xscale-like parameterisation.
  '''

  id_ = 'xscale'

  def __init__(self, params, experiment, reflection, scaled_id=0):
    super(XscaleScaler, self).__init__(params, experiment, reflection, scaled_id)
    (self.g_absorption, self.g_modulation, self.g_decay) = (None, None, None)
    self._initialise_scale_factors()
    self._select_reflections_for_scaling()
    logger.info('Completed initialisation of Xscale Scaler. \n' + '*'*40 + '\n')

  def _initialise_scale_factors(self):
    logger.info('Initialising scale factor objects. \n')
    configdict = self.experiments.scaling_model.configdict
    if 'decay' in configdict['corrections']:
      self._initialise_decay_term()
    if 'absorption' in configdict['corrections']:
      self._initialise_absorption_term()
    if 'modulation' in configdict['corrections']:
      self._initialise_modulation_term()

  def _initialise_decay_term(self):
    '''calculate the 'normalised time', and initialise a SmoothBScaleFactor'''
    refl_table = self.reflection_table
    configdict = self.experiments.scaling_model.configdict
    refl_table['normalised_res_values'] = (((1.0 / (refl_table['d']**2))
      - configdict['resmin']) / configdict['res_bin_width'])
    refl_table['norm_time_values'] = ((refl_table['xyzobs.px.value'].parts()[2]
      - configdict['zmin']) / configdict['time_bin_width'])
    self.g_decay = self.experiments.scaling_model.components['decay']
    self._g_parameterisation['g_decay'] = self.g_decay
    logger.info('The decay term ScaleFactor object was successfully initialised')

  def _initialise_absorption_term(self):
    refl_table = self.reflection_table
    configdict = self.experiments.scaling_model.configdict
    refl_table['normalised_x_abs_values'] = ((refl_table['xyzobs.px.value'].parts()[0]
      - configdict['xmin']) / configdict['x_bin_width'])
    refl_table['normalised_y_abs_values'] = ((refl_table['xyzobs.px.value'].parts()[1]
      - configdict['ymin']) / configdict['y_bin_width'])
    if 'norm_time_values' not in refl_table:
      refl_table['norm_time_values'] = ((refl_table['xyzobs.px.value'].parts()[2]
        - configdict['zmin']) / configdict['time_bin_width'])
    self.g_absorption = self.experiments.scaling_model.components['absorption']
    self._g_parameterisation['g_absorption'] = self.g_absorption
    logger.info('The absorption term ScaleFactor object was successfully initialised')

  def _initialise_modulation_term(self):
    refl_table = self.reflection_table
    configdict = self.experiments.scaling_model.configdict
    refl_table['normalised_x_det_values'] = ((refl_table['xyzobs.px.value'].parts()[0]
      - configdict['xmin']) / configdict['x_det_bin_width'])
    refl_table['normalised_y_det_values'] = ((refl_table['xyzobs.px.value'].parts()[1]
      - configdict['ymin']) / configdict['y_det_bin_width'])
    self.g_modulation = self.experiments.scaling_model.components['modulation']
    self._g_parameterisation['g_modulation'] = self.g_modulation
    logger.info('The modulation term ScaleFactor object was successfully initialised')

  def _select_reflections_for_scaling(self):
    (refl_for_scaling, weights_for_scaling, _) = (
      self._scaling_subset(self.reflection_table, self.params,
      error_model_params=self.params.weighting.error_model_params))
    self._Ih_table = SingleIhTable(refl_for_scaling, weights_for_scaling.weights)
    if self.params.weighting.tukey_biweighting:
      self.Ih_table.apply_tukey_biweighting()
    '''set the normalised values within the scale_factor object to current'''
    if 'modulation' in self.corrections:
      self.g_modulation.update_reflection_data(
        refl_for_scaling['normalised_x_det_values'],
        refl_for_scaling['normalised_y_det_values'])
    if 'decay' in self.corrections:
      self.g_decay.update_reflection_data(
        refl_for_scaling['normalised_res_values'],
        refl_for_scaling['norm_time_values'])
    if 'absorption' in self.experiments.scaling_model.configdict['corrections']:
      self.g_absorption.update_reflection_data(
        refl_for_scaling['normalised_x_abs_values'],
        refl_for_scaling['normalised_y_abs_values'],
        refl_for_scaling['norm_time_values'])

  def apply_selection_to_SFs(self, sel):
    if 'decay' in self.corrections:
      self.g_decay.update_reflection_data(
        normalised_x_values=self.g_decay.normalised_x_values.select(sel),
        normalised_y_values=self.g_decay.normalised_y_values.select(sel))
    if 'absorption' in self.corrections:
      self.g_absorption.update_reflection_data(
        normalised_x_values=self.g_absorption.normalised_x_values.select(sel),
        normalised_y_values=self.g_absorption.normalised_y_values.select(sel),
        normalised_z_values=self.g_absorption.normalised_z_values.select(sel))
    if 'modulation' in self.corrections:
      self.g_modulation.update_reflection_data(
        normalised_x_values=self.g_modulation.normalised_x_values.select(sel),
        normalised_y_values=self.g_modulation.normalised_y_values.select(sel))

  def normalise_scale_component(self):
    pass
    #assert 0, 'method not yet implemented'

  def normalise_decay_component(self):
    pass
    #assert 0, 'method not yet implemented'

  def calc_expanded_scales(self):
    expanded_scale_factors = flex.double([1.0]*len(self.reflection_table))
    if 'modulation' in self.corrections:
      self.g_modulation.update_reflection_data(
        self.reflection_table['normalised_x_det_values'],
        self.reflection_table['normalised_y_det_values'])
      self.g_modulation.calculate_scales()
      expanded_scale_factors *= self.g_modulation.inverse_scales
    if 'decay' in self.corrections:
      self.g_decay.update_reflection_data(
        self.reflection_table['normalised_res_values'],
        self.reflection_table['norm_time_values'])
      self.g_decay.calculate_scales()
      expanded_scale_factors *= self.g_decay.inverse_scales
    if 'absorption' in self.corrections:
      self.g_absorption.update_reflection_data(
        self.reflection_table['normalised_x_abs_values'],
        self.reflection_table['normalised_y_abs_values'],
        self.reflection_table['norm_time_values'])
      self.g_absorption.calculate_scales()
      expanded_scale_factors *= self.g_absorption.inverse_scales
    return expanded_scale_factors


class MultiScaler(ScalerUtilities):
  '''
  Scaler for multiple datasets - takes in params, experiments and
  a list of SingleScalers.
  '''
  def __init__(self, params, experiments, single_scalers):
    '''initialise from a list if single scaler'''
    logger.info('\nInitialising a MultiScaler instance. \n')
    super(MultiScaler, self).__init__()
    self.single_scalers = single_scalers
    self._initial_keys = self.single_scalers[0].initial_keys
    self._params = params
    self._experiments = experiments[0]
    self._Ih_table = JointIhTable(self.single_scalers)
    logger.info('Completed initialisation of MultiScaler. \n' + '*'*40 + '\n')

  def calc_absorption_constraint(self, apm):
    'method only called in aimless scaling'
    R = flex.double([])
    G = flex.double([])
    scaler_ids = [scaler.id_ for scaler in self.single_scalers]
    if 'aimless' in scaler_ids:
      for i, scaler in enumerate(self.single_scalers):
        R.extend(scaler.calc_absorption_constraint(apm.apm_list[i])[0])
        G.extend(scaler.calc_absorption_constraint(apm.apm_list[i])[1])
    return (R, G)

  def update_error_model(self):
    for scaler in self.single_scalers:
      scaler.update_error_model()
    self._Ih_table = JointIhTable(self.single_scalers)

  def update_for_minimisation(self, apm):
    '''update the scale factors and Ih for the next iteration of minimisation,
    update the x values from the amp to the individual apms, as this is where
    basis functions, target functions etc get access to the parameters.'''
    for i, _ in enumerate(apm.n_params_in_each_apm):
      apm.apm_list[i].x = apm.x[apm.n_cumul_params_list[i]:apm.n_cumul_params_list[i+1]]
    apm.active_derivatives = sparse.matrix(self.Ih_table.size, apm.n_active_params)
    for i, scaler in enumerate(self.single_scalers):
      basis_fn = scaler.get_basis_function(apm.apm_list[i])
      scaler.Ih_table.inverse_scale_factors = basis_fn[0]
      expanded = basis_fn[1].transpose() * self.Ih_table.h_index_expand_list[i]
      apm.active_derivatives.assign_block(expanded.transpose(), 0, apm.n_cumul_params_list[i])
    self.Ih_table.calc_Ih()

  def expand_scales_to_all_reflections(self):
    for scaler in self.single_scalers:
      scaler.expand_scales_to_all_reflections(caller=self)

  def join_multiple_datasets(self):
    '''method to create a joint reflection table'''
    joined_reflections = flex.reflection_table()
    joined_outliers = flex.reflection_table()
    for scaler in self.single_scalers:
      joined_reflections.extend(scaler.reflection_table)
      joined_outliers.extend(scaler.outlier_table)
    self._outlier_table = joined_outliers
    miller_set = miller.set(crystal.symmetry(
      space_group=self.single_scalers[0].experiments.crystal.get_space_group()),
      indices=joined_reflections['asu_miller_index'], anomalous_flag=False)
    permuted = miller_set.sort_permutation(by_value='packed_indices')
    self._reflection_table = joined_reflections.select(permuted)
    weights = Weighting(self.reflection_table).weights
    self._Ih_table = SingleIhTable(self.reflection_table, weights)
    self._reflection_table['Ih_values'] = self.Ih_table.Ih_values
    if self.params.scaling_options.reject_outliers:
      sel = reject_outliers(self, self.params.scaling_options.outlier_zmax)
      outlier_sel = ~sel
      self._outlier_table.extend(self._reflection_table.select(outlier_sel))
      n_outliers = sel.count(False)
      msg = ('Combined outlier rejection has been performed across all datasets, {sep}'
        '{0} outliers were found which have been removed from the dataset. {sep}'.format(
        n_outliers, sep='\n'))
      logger.info(msg)
      self._reflection_table = self._reflection_table.select(sel)
      self.Ih_table = self.Ih_table.select(sel)
      self._reflection_table['Ih_values'] = self.Ih_table.Ih_values
      msg = ('A new best estimate for I_h for all reflections across all datasets {sep}'
        'has now been calculated. {sep}').format(sep='\n')
      logger.info(msg)

  def calc_merging_statistics(self):
    joint_result, _ = super(MultiScaler, self).calc_merging_statistics()
    results = []
    scaled_ids = []
    for scaler in self.single_scalers:
      result, scaled_id = scaler.calc_merging_statistics()
      results.append(result[0])
      scaled_ids.append(scaled_id)
    results.append(joint_result[0])
    scaled_ids.append('x')
    return (results, scaled_ids)

class TargetScaler(MultiScaler):
  '''
  Target Scaler for scaling one dataset against already scaled data - takes in
  params, lists of scaled and unscaled experiments, a list of already scaled
  SingleScalers and a list of unscaled reflections.
  '''
  def __init__(self, params, scaled_experiments, scaled_scalers, unscaled_experiments, unscaled_scalers):
    logger.info('\nInitialising a TargetScaler instance. \n')
    super(TargetScaler, self).__init__(params, scaled_experiments, scaled_scalers)
    self.unscaled_scalers = unscaled_scalers
    self._experiments = unscaled_experiments[0]
    #replace above with ScalerFactory to allow for scaling multiple against multiple?
    target_Ih_table = self.Ih_table
    for scaler in unscaled_scalers:
      for i, miller_idx in enumerate(scaler.Ih_table.asu_miller_index):
        sel = target_Ih_table.asu_miller_index == miller_idx
        Ih_values = target_Ih_table.Ih_values.select(sel)
        if Ih_values:
          scaler.Ih_table.Ih_values[i] = Ih_values[0] #all same so just get [0]
        else:
          scaler.Ih_table.Ih_values[i] = 0.0 #set to zero to allow selection below
      sel = scaler.Ih_table.Ih_values != 0.0
      scaler.Ih_table = scaler.Ih_table.select(sel)
      scaler.apply_selection_to_SFs(sel)
    if len(self.unscaled_scalers) > 1:
      self._Ih_table = JointIhTable(self.unscaled_scalers)
    else:
      self._Ih_table = self.unscaled_scalers[0].Ih_table
    self._g_parameterisation = self.unscaled_scalers[0].g_parameterisation
    logger.info('Completed initialisation of TargetScaler. \n' + '*'*40 + '\n')

  def get_target_function(self, apm):
    '''override the target function method for fixed Ih'''
    return target_function_fixedIh(self, apm).return_targets()

  def update_for_minimisation(self, apm):
    '''update the scale factors and Ih for the next iteration of minimisation'''
    if len(self.unscaled_scalers) == 1:
      basis_fn = self.unscaled_scalers[0].get_basis_function(apm)
      apm.active_derivatives = basis_fn[1]
      self.unscaled_scalers[0].Ih_table.inverse_scale_factors = basis_fn[0]
    else:
      '''update the scale factors and Ih for the next iteration of minimisation,
      update the x values from the amp to the individual apms, as this is where
      basis functions, target functions etc get access to the parameters.'''
      for i, _ in enumerate(apm.n_params_in_each_apm):
        apm.apm_list[i].x = apm.x[apm.n_cumul_params_list[i]:apm.n_cumul_params_list[i+1]]
      apm.active_derivatives = sparse.matrix(self.Ih_table.size, apm.n_active_params)
      for i, scaler in enumerate(self.unscaled_scalers):
        basis_fn = scaler.get_basis_function(apm.apm_list[i])
        scaler.Ih_table.inverse_scale_factors = basis_fn[0]
        expanded = basis_fn[1].transpose() * self.Ih_table.h_index_expand_list[i]
        apm.active_derivatives.assign_block(expanded.transpose(), 0, apm.n_cumul_params_list[i])
    #note - we don't calculate Ih here as using a target instead

  def calc_absorption_constraint(self, apm):
    'method only called in aimless scaling'
    R = flex.double([])
    G = flex.double([])
    scaler_ids = [scaler.id_ for scaler in self.unscaled_scalers]
    if 'aimless' in scaler_ids:
      if len(self.unscaled_scalers) == 1:
        R.extend(self.unscaled_scalers[0].calc_absorption_constraint(apm)[0])
        G.extend(self.unscaled_scalers[0].calc_absorption_constraint(apm)[1])
      else:
        for i, scaler in enumerate(self.unscaled_scalers):
          R.extend(scaler.calc_absorption_constraint(apm.apm_list[i])[0])
          G.extend(scaler.calc_absorption_constraint(apm.apm_list[i])[1])
    return (R, G)

  def expand_scales_to_all_reflections(self):
    for scaler in self.unscaled_scalers:
      scaler.expand_scales_to_all_reflections(caller=self)

  def calc_merging_statistics(self):
    return self.dm1.calc_merging_statistics()

  def join_multiple_datasets(self):
    '''method to create a joint reflection table'''
    joined_reflections = flex.reflection_table()
    for scaler in self.single_scalers:
      joined_reflections.extend(scaler.reflection_table)
    for scaler in self.unscaled_scalers:
      joined_reflections.extend(scaler.reflection_table)
    miller_set = miller.set(crystal.symmetry(
      space_group=self.single_scalers[0].experiments.crystal.get_space_group()),
      indices=joined_reflections['asu_miller_index'], anomalous_flag=False)
    permuted = miller_set.sort_permutation(by_value='packed_indices')
    self._reflection_table = joined_reflections.select(permuted)
    #weights = Weighting(self.reflection_table).weights
    #self._Ih_table = SingleIhTable(self.reflection_table, weights)
    #self._reflection_table['Ih_values'] = self.Ih_table.Ih_values
    #if self.params.scaling_options.reject_outliers:
    #  sel = reject_outliers(self, self.params.scaling_options.outlier_zmax)
    #  n_outliers = sel.count(False)
    #  msg = ('Combined outlier rejection has been performed across all datasets, {sep}'
    #    '{0} outliers were found which have been removed from the dataset. {sep}'.format(
    #    n_outliers, sep='\n'))
    #  logger.info(msg)
    #  self._reflection_table = self._reflection_table.select(sel)
    #  self.Ih_table = self.Ih_table.select(sel)
    #  self._reflection_table['Ih_values'] = self.Ih_table.Ih_values
    #  msg = ('A new best estimate for I_h for all reflections across all datasets {sep}'
    #    'has now been calculated. {sep}').format(sep='\n')
    #  logger.info(msg)
