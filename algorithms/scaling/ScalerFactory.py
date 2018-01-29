'''
Collection of factories for creating the scalers.
'''
import logging
import cPickle as pickle
import numpy as np
from dials.array_family import flex
from cctbx import miller, crystal
from scitbx import sparse
import iotbx.merging_statistics
from target_function import (target_function,
  target_function_fixedIh, xds_target_function_log)
import basis_functions as bf
from scaling_utilities import calc_s2d, sph_harm_table
from Wilson_outlier_test import (
  calculate_wilson_outliers, calc_normE2)
from reflection_weighting import Weighting
from target_Ih import SingleIhTable, JointIhTable, IhTableBase
import minimiser_functions as mf
from aimless_outlier_rejection import reject_outliers

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
    single_scalers = targetscaler.data_managers
    single_scalers.append(targetscaler.dm1)
    return MultiScaler(targetscaler.params, [targetscaler.experiments], single_scalers)

class TargetScalerFactory(object):
  'Factory for creating a targeted scaler for multiple datasets'
  @classmethod
  def create(cls, params, experiments, reflections, is_scaled_list):
    '''sort scaled and unscaled datasets to pass to TargetScaler'''
    scaled_experiments = []
    scaled_scalers = []
    unscaled_scalers = []
    for i, reflection in enumerate(reflections):
      if is_scaled_list[i] is True:
        scaled_experiments.append(experiments[i])
        scaled_scalers.append(SingleScalerFactory.create(params, experiments[i],
          reflection, scaled_id=i))
      else:
        unscaled_scalers.append(SingleScalerFactory.create(params, experiments[i],
          reflection, scaled_id=i))
    return TargetScaler(params, scaled_experiments, scaled_scalers, unscaled_scalers[0])


class ScalerUtilities(object):
  '''Base class for all Scalers (single and multiple)'''
  def __init__(self):
    'General attributes relevant for all parameterisations'
    self._experiments = None
    self._params = None
    self._reflection_table = None
    self._Ih_table = None
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
  def params(self):
    return self._params

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
    return bf.basis_function(self, apm).return_basis()

  def expand_scales_to_all_reflections(self, caller=None):
    '''method to be filled in by subclasses'''
    pass

  def clean_reflection_table(self):
    '''method to be filled in by subclasses'''
    pass

  def save_reflection_table(self, filename):
    ''' Save the reflections to file. '''
    self.reflection_table.as_pickle(filename)

  def save_scaler(self, filename):
    ''' Save the data manager to file. '''
    data_file = open(filename, 'w')
    pickle.dump(self, data_file)
    data_file.close()

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
      indices=self.Ih_table.asu_miller_index, anomalous_flag=False)
    scaled_intensities = self.Ih_table.intensities/self.Ih_table.inverse_scale_factors
    sigmas = (self.Ih_table.variances**0.5)/self.Ih_table.inverse_scale_factors
    i_obs = miller.array(miller_set, data=scaled_intensities, sigmas=sigmas)
    i_obs.set_observation_type_xray_intensity()
    result = iotbx.merging_statistics.dataset_statistics(
      i_obs=i_obs, n_bins=20, anomalous=False, use_internal_variance=False,
      eliminate_sys_absent=False, sigma_filtering=None)
    return [result]


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
    self._params = params
    logger.info("Dataset id for this reflection table is %s." % scaled_id)
    logger.info(('The type of scaling model being applied to this dataset {sep}'
      'is an {0}. {sep}').format(self.experiments.scaling_model.id_, sep='\n'))
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
  def initial_keys(self):
    '''list of initial reflection table keys.'''
    return self._initial_keys

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
    reflections = reflections.select(reflections['partiality'] > 0.9)
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


class KBScaler(SingleScaler):
  '''
  Scaler for single dataset using simple KB parameterisation.
  '''
  def __init__(self, params, experiment, reflection, scaled_id=0):
    super(KBScaler, self).__init__(params, experiment, reflection, scaled_id)
    (self.g_scale, self.g_decay) = (None, None)
    self._initialise_scale_factors()
    self._select_reflections_for_scaling()
    logger.info('Completed initialisation of KB Scaler. \n' + '*'*40 + '\n')

  def _initialise_scale_factors(self):
    if self.params.parameterisation.scale_term:
      self.g_scale = self.experiments.scaling_model.components['scale']
      self._g_parameterisation['g_scale'] = self.g_scale
    if self.params.parameterisation.decay_term:
      self.g_decay = self.experiments.scaling_model.components['decay']
      self._g_parameterisation['g_decay'] = self.g_decay

  def _select_reflections_for_scaling(self):
    (reflections_for_scaling, weights_for_scaling, _) = (
      self._scaling_subset(self.reflection_table, self.params))
    self._Ih_table = SingleIhTable(reflections_for_scaling, weights_for_scaling.weights)
    if self.params.parameterisation.scale_term:
      self.g_scale.update_reflection_data(n_refl=self.Ih_table.size)
    if self.params.parameterisation.decay_term:
      self.g_decay.update_reflection_data(dvalues=reflections_for_scaling['d'])

  def apply_selection_to_SFs(self, sel):
    '''Updates data from within current SF objects using the given selection.
       Required for targeted scaling.'''
    if self.params.parameterisation.decay_term:
      self.g_decay.update_reflection_data(dvalues=self.g_decay.d_values.select(sel))
    if self.params.parameterisation.scale_term:
      self.g_scale.update_reflection_data(n_refl=sel.count(True))

  def expand_scales_to_all_reflections(self, caller=None):
    expanded_scale_factors = flex.double([1.0]*len(self.reflection_table))
    if self.params.parameterisation.scale_term:
      self.g_scale.update_reflection_data(n_refl=len(self.reflection_table))
      self.g_scale.calculate_scales_and_derivatives()
      expanded_scale_factors *= self.g_scale.inverse_scales
      logger.info(('Scale factor K = {0:.4f} determined during minimisation {sep}'
        'has now been applied to all reflections. {sep}').format(
        list(self.g_scale.parameters)[0], sep='\n'))
    if self.params.parameterisation.decay_term:
      self.g_decay.update_reflection_data(dvalues=self.reflection_table['d'])
      self.g_decay.calculate_scales_and_derivatives()
      expanded_scale_factors *= self.g_decay.inverse_scales
      logger.info(('B-factor B = {0:.4f} determined during minimisation {sep}'
        'has now been applied to all reflections. {sep}').format(
        list(self.g_decay.parameters)[0], sep='\n'))
    self._reflection_table['inverse_scale_factor'] = expanded_scale_factors
    self._Ih_table = SingleIhTable(self.reflection_table)
    weights = self._update_weights_for_scaling(self.reflection_table,
      self.params, weights_filter=False, error_model_params=None).weights
    self.Ih_table.weights = weights
    self.Ih_table = self.Ih_table.select(self.Ih_table.weights != 0.0)
    self._reflection_table = self._reflection_table.select(weights != 0.0)
    self._reflection_table['Ih_values'] = self.Ih_table.Ih_values


class AimlessScaler(SingleScaler):
  '''
  Scaler for single dataset using aimless-like parameterisation.
  '''
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

  def round_of_outlier_rejection(self):
    '''calculate outliers from the reflections in the Ih_table,
    and use these to filter the reflection table and Ih_table.'''
    sel = reject_outliers(self, self.params.scaling_options.outlier_zmax)
    self._reflection_table = self._reflection_table.select(sel)
    self.Ih_table = self.Ih_table.select(sel)
    self._reflection_table['Ih_values'] = self.Ih_table.Ih_values
    msg = ('A round of outlier rejection has been performed, {0} outliers {sep}'
        'were found which have been removed from the dataset. {sep}'.format(
        sel.count(False), sep='\n'))
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
    self.g_scale.update_reflection_data(refl_for_scaling['norm_rot_angle'])
    if self.params.parameterisation.decay_term:
      self.g_decay.update_reflection_data(dvalues=refl_for_scaling['d'],
        normalised_values=refl_for_scaling['norm_time_values'])
    if self.params.parameterisation.absorption_term:
      sph_harm_table_T = self.sph_harm_table.transpose()
      sel_sph_harm_table = sph_harm_table_T.select_columns(selection.iselection())
      self.g_absorption.update_reflection_data(sel_sph_harm_table.transpose())

  def apply_selection_to_SFs(self, sel):
    '''Updates data from within current SF objects using the given selection.
       Required for targeted scaling.'''
    self.g_scale.update_reflection_data(self.g_scale.normalised_values.select(sel))
    if self.params.parameterisation.decay_term:
      self.g_decay.update_reflection_data(dvalues=self.g_decay.d_values.select(sel),
        normalised_values=self.g_decay.normalised_values.select(sel))
    if self.params.parameterisation.absorption_term:
      sph_harm_table_T = self.g_absorption.harmonic_values.transpose()
      sel_sph_harm_table = sph_harm_table_T.select_columns(sel.iselection())
      self.g_absorption.update_reflection_data(sel_sph_harm_table.transpose())

  def update_error_model(self):
    '''apply a correction to try to improve the error estimate.'''
    self.params.weighting.error_model_params = (
      mf.error_scale_LBFGSoptimiser(self.Ih_table, flex.double([1.0, 0.01])).x)
    self.Ih_table.update_aimless_error_model(self.params.weighting.error_model_params)

  def _initialise_scale_factors(self):
    '''initialise scale factors and add to self.active_parameters'''
    #logger.info('Initialising scale factor objects. \n')
    self._initialise_scale_term(self.reflection_table)
    if self.params.parameterisation.decay_term:
      self._initialise_decay_term(self.reflection_table)
    if self.params.parameterisation.absorption_term:
      self._initialise_absorption_term(self.reflection_table,
        self.params.parameterisation.lmax)
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
    '''if 'scaling' in self.experiments.__dict__:
      if self.experiments.scaling.absorption.object_type == "SHScaleFactor":
        parameters = self.experiments.scaling.absorption.parameters
        self.g_absorption = SF.SHScaleFactor(parameters, n_param)
    else:
      self.g_absorption = SF.SHScaleFactor(0.0, n_param)'''
    self.g_absorption = self.experiments.scaling_model.components['absorption']
    self._g_parameterisation['g_absorption'] = self.g_absorption
    msg = ('The absorption term ScaleFactor object was successfully initialised. {sep}'
      'The absorption term will be parameterised by a set of spherical {sep}'
      'harmonics up to an lmax of {0} ({1} parameters). {sep}'
      ).format(lmax, self.g_absorption.n_params, sep='\n')
    logger.info(msg)

  def calc_absorption_constraint(self, apm):
    #should only be called from target function if g_absorption in active params
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

  def _normalise_scales_and_B(self):
    if self.params.parameterisation.decay_term:
      maxB = max(flex.double(np.log(self.g_decay.inverse_scales))
                 * 2.0 * (self.g_decay.d_values**2))
      self.g_decay.parameters -= flex.double([maxB] * self.g_decay.n_params)
      self.g_decay.calculate_scales_and_derivatives()
    sel = (self.g_scale.normalised_values == min(self.g_scale.normalised_values))
    initial_scale = self.g_scale.inverse_scales.select(sel)[0]
    self.g_scale.parameters /= initial_scale
    self.g_scale.calculate_scales_and_derivatives()

  def expand_scales_to_all_reflections(self, caller=None):
    expanded_scale_factors = flex.double([1.0]*len(self.reflection_table))
    self.g_scale.update_reflection_data(self.reflection_table['norm_rot_angle'])
    self.g_scale.calculate_scales()
    expanded_scale_factors *= self.g_scale.inverse_scales
    if self.params.parameterisation.decay_term:
      self.g_decay.update_reflection_data(dvalues=self.reflection_table['d'],
        normalised_values=self.reflection_table['norm_time_values'])
      self.g_decay.calculate_scales()
      expanded_scale_factors *= self.g_decay.inverse_scales
    if self.params.parameterisation.absorption_term:
      self.g_absorption.update_reflection_data(self.sph_harm_table)
      self.g_absorption.calculate_scales_and_derivatives()
      expanded_scale_factors *= self.g_absorption.inverse_scales
    self._reflection_table['inverse_scale_factor'] = expanded_scale_factors
    logger.info(('Scale factors determined during minimisation have now been applied {sep}'
      'to all reflections. {sep}').format(sep='\n'))
    self._Ih_table = SingleIhTable(self.reflection_table)
    weights = self._update_weights_for_scaling(self.reflection_table,
      self.params, weights_filter=False, error_model_params=None).weights
    self.Ih_table.weights = weights
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

  def clean_reflection_table(self):
    self._initial_keys.append('inverse_scale_factor')
    self._initial_keys.append('Ih_values')
    #keep phi column for now for comparing to aimless
    self._initial_keys.append('phi')
    for key in self.reflection_table.keys():
      if not key in self._initial_keys:
        del self.reflection_table[key]


class XscaleScaler(SingleScaler):
  '''
  Scaler for single dataset using xscale-like parameterisation.
  '''
  def __init__(self, params, experiment, reflection, scaled_id=0):
    super(XscaleScaler, self).__init__(params, experiment, reflection, scaled_id)
    assert 0, 'method not yet implemented'


class MultiScaler(ScalerUtilities):
  '''
  Scaler for multiple datasets - takes in params, experiments and
  a list of SingleScalers.
  '''
  def __init__(self, params, experiments, single_scalers):
    '''initialise from a list if single scaler'''
    logger.info('\nInitialising a MultiScaler instance. \n')
    super(MultiScaler, self).__init__()
    self.data_managers = single_scalers
    self._params = params
    self._experiments = experiments[0]
    self._Ih_table = JointIhTable(self.data_managers)
    logger.info('Completed initialisation of MultiScaler. \n' + '*'*40 + '\n')

  def calc_absorption_constraint(self, apm):
    'method only called in aimless scaling'
    R = flex.double([])
    G = flex.double([])
    for i, dm in enumerate(self.data_managers):
      R.extend(dm.calc_absorption_constraint(apm.apm_list[i])[0])
      G.extend(dm.calc_absorption_constraint(apm.apm_list[i])[1])
    return (R, G)

  def update_error_model(self):
    for dm in self.data_managers:
      dm.update_error_model()
    self._Ih_table = JointIhTable(self.data_managers)

  def update_for_minimisation(self, apm):
    '''update the scale factors and Ih for the next iteration of minimisation,
    update the x values from the amp to the individual apms, as this is where
    basis functions, target functions etc get access to the parameters.'''
    for i, _ in enumerate(apm.n_params_in_each_apm):
      apm.apm_list[i].x = apm.x[apm.n_cumul_params_list[i]:apm.n_cumul_params_list[i+1]]
    apm.active_derivatives = sparse.matrix(self.Ih_table.size, apm.n_active_params)
    for i, dm in enumerate(self.data_managers):
      basis_fn = dm.get_basis_function(apm.apm_list[i])
      dm.Ih_table.inverse_scale_factors = basis_fn[0]
      expanded = basis_fn[1].transpose() * self.Ih_table.h_index_expand_list[i]
      apm.active_derivatives.assign_block(expanded.transpose(), 0, apm.n_cumul_params_list[i])
    self.Ih_table.calc_Ih()

  def expand_scales_to_all_reflections(self, caller=None):
    for dm in self.data_managers:
      dm.expand_scales_to_all_reflections(caller=self)

  def join_multiple_datasets(self):
    '''method to create a joint reflection table'''
    joined_reflections = flex.reflection_table()
    for dm in self.data_managers:
      joined_reflections.extend(dm.reflection_table)
    miller_set = miller.set(crystal.symmetry(
      space_group=self.data_managers[0].experiments.crystal.get_space_group()),
      indices=joined_reflections['asu_miller_index'], anomalous_flag=True)
    permuted = miller_set.sort_permutation(by_value='packed_indices')
    self._reflection_table = joined_reflections.select(permuted)
    weights = Weighting(self.reflection_table).weights
    self._Ih_table = SingleIhTable(self.reflection_table, weights)
    self._reflection_table['Ih_values'] = self.Ih_table.Ih_values
    if self.params.scaling_options.reject_outliers:
      sel = reject_outliers(self, self.params.scaling_options.outlier_zmax)
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
    joint_result = super(MultiScaler, self).calc_merging_statistics()[0]
    results = []
    for dm in self.data_managers:
      results.append(dm.calc_merging_statistics()[0])
    results.append(joint_result)
    return results

  def clean_reflection_table(self):
    self.data_managers[0].initial_keys.append('inverse_scale_factor')
    self.data_managers[0].initial_keys.append('Ih_values')
    #keep phi column for now for comparing to aimless
    self.data_managers[0].initial_keys.append('phi')
    for key in self.data_managers[0].reflection_table.keys():
      if not key in self.data_managers[0].initial_keys:
        del self.reflection_table[key]


class TargetScaler(MultiScaler):
  '''
  Target Scaler for scaling one dataset against already scaled data - takes in
  params, lists of scaled and unscaled experiments, a list of already scaled
  SingleScalers and a list of unscaled reflections.
  '''
  def __init__(self, params, scaled_experiments, scaled_scalers, unscaled_scaler):
    logger.info('\nInitialising a TargetScaler instance. \n')
    super(TargetScaler, self).__init__(params, scaled_experiments, scaled_scalers)
    self.dm1 = unscaled_scaler
    #replace above with ScalerFactory to allow for scaling multiple against multiple?
    target_Ih_table = self.Ih_table
    for i, miller_idx in enumerate(self.dm1.Ih_table.asu_miller_index):
      sel = target_Ih_table.asu_miller_index == miller_idx
      Ih_values = target_Ih_table.Ih_values.select(sel)
      if Ih_values:
        self.dm1.Ih_table.Ih_values[i] = Ih_values[0] #all same so just get [0]
      else:
        self.dm1.Ih_table.Ih_values[i] = 0.0 #set to zero to allow selection below
    sel = self.dm1.Ih_table.Ih_values != 0.0
    self.dm1.Ih_table = self.dm1.Ih_table.select(sel)
    self.dm1.apply_selection_to_SFs(sel)
    self._g_parameterisation = self.dm1.g_parameterisation
    logger.info('Completed initialisation of TargetScaler. \n' + '*'*40 + '\n')

  def get_target_function(self, apm):
    '''override the target function method for fixed Ih'''
    return target_function_fixedIh(self.dm1, apm).return_targets()

  def update_for_minimisation(self, apm):
    '''update the scale factors and Ih for the next iteration of minimisation'''
    basis_fn = self.dm1.get_basis_function(apm)
    apm.active_derivatives = basis_fn[1]
    self.dm1.Ih_table.inverse_scale_factors = basis_fn[0]
    #note - we don't calculate Ih here as using a target instead

  def expand_scales_to_all_reflections(self, caller=None):
    return self.dm1.expand_scales_to_all_reflections(caller=self)

  def calc_merging_statistics(self):
    return self.dm1.calc_merging_statistics()

  def clean_reflection_table(self):
    self.dm1.clean_reflection_table()

  #def save_reflection_table(self, filename):
  #  ''' Save the reflections to file. '''
  #  self.dm1.save_reflection_table(filename)#='target_scaled.pickle')

  def join_multiple_datasets(self):
    '''method to create a joint reflection table'''
    joined_reflections = flex.reflection_table()
    for scaler in self.data_managers:
      joined_reflections.extend(scaler.reflection_table)
    joined_reflections.extend(self.dm1.reflection_table)
    miller_set = miller.set(crystal.symmetry(
      space_group=self.data_managers[0].experiments.crystal.get_space_group()),
      indices=joined_reflections['asu_miller_index'], anomalous_flag=True)
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
