import abc
import logging
import numpy as np
from dials.array_family import flex
from cctbx import miller, crystal
from scitbx import sparse
from dials_scaling_helpers_ext import row_multiply
import iotbx.merging_statistics
from dials.algorithms.scaling.target_function import \
  target_function, target_function_fixedIh
from dials.algorithms.scaling.basis_functions import basis_function
from dials.algorithms.scaling.scaling_utilities import (sph_harm_table,
  reject_outliers, calculate_wilson_outliers, calc_normE2)
from dials.algorithms.scaling.reflection_weighting import Weighting
from dials.algorithms.scaling.Ih_table import SingleIhTable, JointIhTable, IhTableBase
from dials.algorithms.scaling.minimiser_functions import error_scale_LBFGSoptimiser
logger = logging.getLogger('dials')

class ScalerBase(object):
  '''Base class for all Scalers (single and multiple)'''

  __metaclass__ = abc.ABCMeta

  def __init__(self):
    'General attributes relevant for all parameterisations'
    self._experiments = None
    self._params = None
    self._reflection_table = []
    self._outlier_table = flex.reflection_table()
    self._Ih_table = None
    self._var_cov = None
    self._initial_keys = []

  @property
  def Ih_table(self):
    return self._Ih_table

  @Ih_table.setter
  def Ih_table(self, new_Ih_table):
    assert isinstance(new_Ih_table, IhTableBase)
    self._Ih_table = new_Ih_table

  @property
  def var_cov_matrix(self):
    return self._var_cov

  @var_cov_matrix.setter
  def var_cov_matrix(self, var_cov):
    self._var_cov = var_cov

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
    self._initial_keys.append('inverse_scale_factor_variance')
    self._initial_keys.append('Ih_values')
    for key in self.reflection_table.keys():
      if not key in self._initial_keys:
        del self._reflection_table[key]

  @abc.abstractmethod
  def update_for_minimisation(self, apm):
    '''update the scale factors and Ih for the next iteration of minimisation'''
    pass

  @abc.abstractmethod
  def expand_scales_to_all_reflections(self, caller=None):
    '''expand scales from a subset to all reflections'''
    pass

  def get_target_function(self, apm):
    '''call the target function'''
    return target_function(self, apm).return_targets()

  def get_residuals_jacobian_weight(self, apm):
    return target_function(self, apm).return_residuals_jacobian_weight()

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
      indices=self.reflection_table['miller_index'], anomalous_flag=False)
    scaled_intensities = (self.reflection_table['intensity']/
      self.reflection_table['inverse_scale_factor'])
    sigmas = ((self.reflection_table['variance']**0.5)/
      self.reflection_table['inverse_scale_factor'])
    i_obs = miller.array(miller_set, data=scaled_intensities)
    i_obs.set_observation_type_xray_intensity()
    i_obs.set_sigmas(sigmas)
    scaled_ids = list(set(self.reflection_table['id']))
    result = iotbx.merging_statistics.dataset_statistics(
      i_obs=i_obs, n_bins=20, anomalous=False, sigma_filtering=None,
      use_internal_variance=True, eliminate_sys_absent=False)
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



class SingleScalerBase(ScalerBase):
  '''
  Parent class for single dataset Scalers, containing a standard
  setup routine for the reflection_table - takes in params, experiment
  and reflection.
  '''

  def __init__(self, params, experiment, reflection, scaled_id=0):
    logger.info('\nInitialising a Single Scaler instance. \n')
    super(SingleScalerBase, self).__init__()
    self._experiments = experiment
    self._params = params
    logger.info("Dataset id for this reflection table is %s." % scaled_id)
    logger.info(('The type of scaling model being applied to this dataset {sep}'
      'is {0}. {sep}').format(self.experiments.scaling_model.id_, sep='\n'))
    #rename the id for now so that we have unique dataset ids
    reflection['id'] = flex.int([scaled_id]*len(reflection))
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
  def components(self):
    'shortcut to scaling model components'
    return self.experiments.scaling_model.components

  @abc.abstractproperty
  def consecutive_scaling_order(self):
    '''should return a nested list of correction names, to indicate the order
    to perform scaling in consecutive scaling mode if concurrent_scaling=0.
    e.g. [['scale', 'decay'], ['absorption']] would cause the first cycle to
    refine scale and decay, and then absorption in a subsequent cycle.'''
    pass

  def update_for_minimisation(self, apm):
    '''update the scale factors and Ih for the next iteration of minimisation'''
    basis_fn = self.get_basis_function(apm)
    apm.derivatives = basis_fn[1]
    self.Ih_table.inverse_scale_factors = basis_fn[0]
    self.Ih_table.calc_Ih()

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
    if self.var_cov_matrix:
      self._reflection_table['inverse_scale_factor_variance'] = self.calc_sf_variances()
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
      not isinstance(caller, MultiScalerBase)):
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

  @abc.abstractmethod
  def apply_selection_to_SFs(self, sel):
    '''Updates data from within current SF objects using the given selection.
       Required for targeted scaling. To be filled in by subclasses'''
    pass

  @abc.abstractmethod
  def calc_expanded_scales(self):
    '''calculate the scale factors for all reflections from the model.
    To be filled in by subclasses.'''
    pass

  def normalise_scale_component(self):
    '''can be filled in by subclass to allow normalising, as is the case for aimless'''
    pass

  def normalise_decay_component(self):
    '''can be filled in by subclass to allow normalising, as is the case for aimless'''
    pass

  def print_scale_init_msg(self):
    msg = ('The scale factor corrections have been successfully initialised {sep}'
      'for the following corrections;{sep}'
      'correction name:      {0} {sep}'
      'number of parameters: {1} {sep}'
        ).format(''.join(str(i)+' ' for i in self.components.iterkeys()),
        ''.join(str(val.n_params)+' ' for val in self.components.itervalues()),
        sep='\n')
    logger.info(msg)

  def calc_absorption_constraint(self, apm):
    '''calculates a constraint for the spherical harmonic absorption correction.
    Should only be called from target function if absorption in active params.'''
    if self.id_ == 'aimless':
      if 'absorption' in apm.components:
        abs_params = apm.select_parameters('absorption')
        residual = (self.absorption_weights * (abs_params)**2)
        gradient = (2.0 * self.absorption_weights * abs_params)
        #return a gradient vector to be added to that calculated in target function
        gradient_vector = flex.double([])
        for comp in apm.components:
          if comp != 'absorption':
            gradient_vector.extend(flex.double([0.0] * apm.components[comp]['n_params']))
          elif comp == 'absorption':
            gradient_vector.extend(gradient)
        return (residual, gradient_vector)
      else:
        return (flex.double([0.0]), flex.double([0.0] * apm.n_active_params))
    return (flex.double([0.0]), flex.double([0.0] * apm.n_active_params))

  def calc_sf_variances(self):
    '''use the parameter var_cov matrix to calculate the variances of the inverse scales'''
    n_param = 0
    for component in self.components:
      n_param += self.components[component].n_params
      n_refl = len(self.components[component].inverse_scales) #should all be same
    jacobian = sparse.matrix(n_refl, n_param)
    n_cumulative_param = 0
    for component in self.components:
      block = self.components[component].derivatives
      n_param = self.components[component].n_params
      for component_2 in self.components:
        if component_2 != component:
          block = row_multiply(block, self.components[component].inverse_scales)
      jacobian.assign_block(block, 0, n_cumulative_param)
      n_cumulative_param += n_param
    jacobian_transpose = jacobian.transpose()

    length = int(len(self.var_cov_matrix)**0.5)
    var_cov_mat = sparse.matrix(length, length)
    for i in range(length):
      col = sparse.matrix_column(length)
      for j in range(length):
        col[j] = self.var_cov_matrix[(i*length) + j]
      var_cov_mat[:, i] = col
    logger.info('Calculating error estimates of inverse scale factors. \n')
    sigmasq = flex.float([])
    #note: must be a faster way to do this next bit?
    for col in jacobian_transpose.cols(): #iterating over reflections
      a = flex.double(col.as_dense_vector())
      var = (a * var_cov_mat) * a
      sigmasq.append(flex.sum(var))
    return sigmasq.as_double()

class KBScaler(SingleScalerBase):
  '''
  Scaler for single dataset using simple KB parameterisation.
  '''

  id_ = 'KB'

  def __init__(self, params, experiment, reflection, scaled_id=0):
    super(KBScaler, self).__init__(params, experiment, reflection, scaled_id)
    self.print_scale_init_msg()
    self._select_reflections_for_scaling()
    logger.info('Completed initialisation of KB Scaler. \n' + '*'*40 + '\n')

  @property
  def consecutive_scaling_order(self):
    return [['scale', 'decay']]

  def _select_reflections_for_scaling(self):
    (refl_for_scaling, weights_for_scaling, _) = (
      self._scaling_subset(self.reflection_table, self.params))
    self._Ih_table = SingleIhTable(refl_for_scaling, weights_for_scaling.weights)
    if 'scale' in self.components:
      self.components['scale'].update_reflection_data(n_refl=self.Ih_table.size)
    if 'decay' in self.components:
      self.components['decay'].update_reflection_data(dvalues=refl_for_scaling['d'])

  def apply_selection_to_SFs(self, sel):
    '''Updates data from within current SF objects using the given selection.
       Required for targeted scaling.'''
    if 'decay' in self.components:
      self.components['decay'].update_reflection_data(
        dvalues=self.components['decay'].d_values.select(sel))
    if 'scale' in self.components:
      self.components['scale'].update_reflection_data(n_refl=sel.count(True))

  def calc_expanded_scales(self):
    expanded_scale_factors = flex.double([1.0]*len(self.reflection_table))
    if 'scale' in self.components:
      self.components['scale'].update_reflection_data(n_refl=len(self.reflection_table))
      self.components['scale'].calculate_scales_and_derivatives()
      expanded_scale_factors *= self.components['scale'].inverse_scales
      logger.info(('Scale factor K = {0:.4f} determined during minimisation {sep}'
        'has now been applied to all reflections. {sep}').format(
        list(self.components['scale'].parameters)[0], sep='\n'))
    if 'decay' in self.components:
      self.components['decay'].update_reflection_data(dvalues=self.reflection_table['d'])
      self.components['decay'].calculate_scales_and_derivatives()
      expanded_scale_factors *= self.components['decay'].inverse_scales
      logger.info(('B-factor B = {0:.4f} determined during minimisation {sep}'
        'has now been applied to all reflections. {sep}').format(
        list(self.components['decay'].parameters)[0], sep='\n'))
    return expanded_scale_factors

  


class AimlessScaler(SingleScalerBase):
  '''
  Scaler for single dataset using aimless-like parameterisation.
  '''

  id_ = 'aimless'

  def __init__(self, params, experiment, reflection, scaled_id=0):
    super(AimlessScaler, self).__init__(params, experiment, reflection, scaled_id)
    self.sph_harm_table = None
    self.absorption_weights = None
    if self.params.scaling_options.reject_outliers:
      self._Ih_table = SingleIhTable(self.reflection_table)
      self.round_of_outlier_rejection()
    self._initialise_scale_factors()
    self._select_reflections_for_scaling()
    logger.info('Completed initialisation of AimlessScaler. \n' + '*'*40 + '\n')

  @property
  def consecutive_scaling_order(self):
    return [['scale', 'decay'], ['absorption']]

  def _initialise_scale_factors(self):
    '''initialise scale factors and add to self.active_parameters'''
    refl_table = self.reflection_table
    if 'scale' in self.components:
      refl_table['norm_rot_angle'] = (refl_table['xyzobs.px.value'].parts()[2]
        * self.experiments.scaling_model.scale_normalisation_factor)
    if 'decay' in self.components:
      refl_table['norm_time_values'] = (refl_table['xyzobs.px.value'].parts()[2]
        * self.experiments.scaling_model.decay_normalisation_factor)
    refl_table['phi'] = (refl_table['xyzobs.px.value'].parts()[2]
      * self.experiments.scan.get_oscillation()[1])
    if 'absorption' in self.components:
      lmax = self.experiments.scaling_model.configdict['lmax']
      self.sph_harm_table = sph_harm_table(refl_table, self.experiments, lmax)
      self.absorption_weights = flex.double([])
      for i in range(1, lmax+1):
        self.absorption_weights.extend(flex.double([1.0] * ((2*i)+1)))
      self.absorption_weights *= self.params.parameterisation.surface_weight
    self.print_scale_init_msg()

  def _select_reflections_for_scaling(self):
    (refl_for_scaling, weights_for_scaling, selection) = (
      self._scaling_subset(self.reflection_table, self.params,
      error_model_params=self.params.weighting.error_model_params))
    self._Ih_table = SingleIhTable(refl_for_scaling, weights_for_scaling.weights)
    if self.params.weighting.tukey_biweighting:
      self.Ih_table.apply_tukey_biweighting()
    '''refactor the next two operations into extract_reflections?
    reset the normalised values within the scale_factor object to current'''
    if 'scale' in self.components:
      self.components['scale'].update_reflection_data(refl_for_scaling['norm_rot_angle'])
    if 'decay' in self.components:
      self.components['decay'].update_reflection_data(dvalues=refl_for_scaling['d'],
        normalised_values=refl_for_scaling['norm_time_values'])
    if 'absorption' in self.components:
      sph_harm_table_T = self.sph_harm_table.transpose()
      sel_sph_harm_table = sph_harm_table_T.select_columns(selection.iselection())
      self.components['absorption'].update_reflection_data(sel_sph_harm_table.transpose())

  def apply_selection_to_SFs(self, sel):
    '''Updates data from within current SF objects using the given selection.
       Required for targeted scaling.'''
    if 'scale' in self.components:
      self.components['scale'].update_reflection_data(
        self.components['scale'].normalised_values.select(sel))
    if 'decay' in self.components:
      self.components['decay'].update_reflection_data(
        dvalues=self.components['decay'].d_values.select(sel),
        normalised_values=self.components['decay'].normalised_values.select(sel))
    if 'absorption' in self.components:
      sph_harm_table_T = self.components['absorption'].harmonic_values.transpose()
      sel_sph_harm_table = sph_harm_table_T.select_columns(sel.iselection())
      self.components['absorption'].update_reflection_data(sel_sph_harm_table.transpose())

  def normalise_scale_component(self):
    '''Method to do an invariant rescale of the scale at t=0 to one.'''
    sel = (self.components['scale'].normalised_values ==
      min(self.components['scale'].normalised_values))
    initial_scale = self.components['scale'].inverse_scales.select(sel)[0]
    self.components['scale'].parameters /= initial_scale
    self.components['scale'].calculate_scales_and_derivatives()
    logger.info('Rescaled the scale component so that the initial scale is 1.\n')

  def normalise_decay_component(self):
    '''Method to do an invariant rescale of the max B to zero.'''
    maxB = max(flex.double(np.log(self.components['decay'].inverse_scales))
                 * 2.0 * (self.components['decay'].d_values**2))
    self.components['decay'].parameters -= flex.double([maxB] * self.components['decay'].n_params)
    self.components['decay'].calculate_scales_and_derivatives()
    logger.info('Rescaled the decay component so that the max B is 0.\n')

  def calc_expanded_scales(self):
    expanded_scale_factors = flex.double([1.0]*len(self.reflection_table))
    if 'scale' in self.components:
      self.components['scale'].update_reflection_data(self.reflection_table['norm_rot_angle'])
      self.components['scale'].calculate_scales_and_derivatives()
      expanded_scale_factors *= self.components['scale'].inverse_scales
    if 'decay' in self.components:
      self.components['decay'].update_reflection_data(dvalues=self.reflection_table['d'],
        normalised_values=self.reflection_table['norm_time_values'])
      self.components['decay'].calculate_scales_and_derivatives()
      expanded_scale_factors *= self.components['decay'].inverse_scales
    if 'absorption' in self.components:
      self.components['absorption'].update_reflection_data(self.sph_harm_table)
      self.components['absorption'].calculate_scales_and_derivatives()
      expanded_scale_factors *= self.components['absorption'].inverse_scales
    return expanded_scale_factors

  def clean_reflection_table(self):
    self._initial_keys.append('phi')
    super(AimlessScaler, self).clean_reflection_table()


class XscaleScaler(SingleScalerBase):
  '''
  Scaler for single dataset using xscale-like parameterisation.
  '''

  id_ = 'xscale'

  def __init__(self, params, experiment, reflection, scaled_id=0):
    super(XscaleScaler, self).__init__(params, experiment, reflection, scaled_id)
    self._initialise_scale_factors()
    self._select_reflections_for_scaling()
    logger.info('Completed initialisation of Xscale Scaler. \n' + '*'*40 + '\n')

  @property
  def consecutive_scaling_order(self):
    return [['decay'], ['absorption'], ['modulation']]

  def _initialise_scale_factors(self):
    logger.info('Initialising scale factor objects. \n')
    if 'decay' in self.components:
      self._initialise_decay_term()
    if 'absorption' in self.components:
      self._initialise_absorption_term()
    if 'modulation' in self.components:
      self._initialise_modulation_term()
    self.print_scale_init_msg()

  def _initialise_decay_term(self):
    '''calculate the 'normalised time', and initialise a SmoothBScaleFactor'''
    refl_table = self.reflection_table
    configdict = self.experiments.scaling_model.configdict
    refl_table['normalised_res_values'] = (((1.0 / (refl_table['d']**2))
      - configdict['resmin']) / configdict['res_bin_width'])
    refl_table['norm_time_values'] = ((refl_table['xyzobs.px.value'].parts()[2]
      - configdict['zmin']) / configdict['time_bin_width'])

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

  def _initialise_modulation_term(self):
    refl_table = self.reflection_table
    configdict = self.experiments.scaling_model.configdict
    refl_table['normalised_x_det_values'] = ((refl_table['xyzobs.px.value'].parts()[0]
      - configdict['xmin']) / configdict['x_det_bin_width'])
    refl_table['normalised_y_det_values'] = ((refl_table['xyzobs.px.value'].parts()[1]
      - configdict['ymin']) / configdict['y_det_bin_width'])

  def _select_reflections_for_scaling(self):
    (refl_for_scaling, weights_for_scaling, _) = (
      self._scaling_subset(self.reflection_table, self.params,
      error_model_params=self.params.weighting.error_model_params))
    self._Ih_table = SingleIhTable(refl_for_scaling, weights_for_scaling.weights)
    if self.params.weighting.tukey_biweighting:
      self.Ih_table.apply_tukey_biweighting()
    '''set the normalised values within the scale_factor object to current'''
    if 'modulation' in self.components:
      self.components['modulation'].update_reflection_data(
        refl_for_scaling['normalised_x_det_values'],
        refl_for_scaling['normalised_y_det_values'])
    if 'decay' in self.components:
      self.components['decay'].update_reflection_data(
        refl_for_scaling['normalised_res_values'],
        refl_for_scaling['norm_time_values'])
    if 'absorption' in self.components:
      self.components['absorption'].update_reflection_data(
        refl_for_scaling['normalised_x_abs_values'],
        refl_for_scaling['normalised_y_abs_values'],
        refl_for_scaling['norm_time_values'])

  def apply_selection_to_SFs(self, sel):
    if 'decay' in self.components:
      self.components['decay'].update_reflection_data(
        normalised_x_values=self.components['decay'].normalised_x_values.select(sel),
        normalised_y_values=self.components['decay'].normalised_y_values.select(sel))
    if 'absorption' in self.components:
      self.components['absorption'].update_reflection_data(
        normalised_x_values=self.components['absorption'].normalised_x_values.select(sel),
        normalised_y_values=self.components['absorption'].normalised_y_values.select(sel),
        normalised_z_values=self.components['absorption'].normalised_z_values.select(sel))
    if 'modulation' in self.components:
      self.components['modulation'].update_reflection_data(
        normalised_x_values=self.components['modulation'].normalised_x_values.select(sel),
        normalised_y_values=self.components['modulation'].normalised_y_values.select(sel))

  def calc_expanded_scales(self):
    expanded_scale_factors = flex.double([1.0]*len(self.reflection_table))
    if 'modulation' in self.components:
      self.components['modulation'].update_reflection_data(
        self.reflection_table['normalised_x_det_values'],
        self.reflection_table['normalised_y_det_values'])
      self.components['modulation'].calculate_scales_and_derivatives()
      expanded_scale_factors *= self.components['modulation'].inverse_scales
    if 'decay' in self.components:
      self.components['decay'].update_reflection_data(
        self.reflection_table['normalised_res_values'],
        self.reflection_table['norm_time_values'])
      self.components['decay'].calculate_scales_and_derivatives()
      expanded_scale_factors *= self.components['decay'].inverse_scales
    if 'absorption' in self.components:
      self.components['absorption'].update_reflection_data(
        self.reflection_table['normalised_x_abs_values'],
        self.reflection_table['normalised_y_abs_values'],
        self.reflection_table['norm_time_values'])
      self.components['absorption'].calculate_scales_and_derivatives()
      expanded_scale_factors *= self.components['absorption'].inverse_scales
    return expanded_scale_factors

class MultiScalerBase(ScalerBase):
  '''Base class for Scalers handling multiple datasets'''
  def __init__(self, params, experiments, single_scalers):
    '''initialise from a list of single scalers'''
    super(MultiScalerBase, self).__init__()
    self.single_scalers = single_scalers
    self._initial_keys = self.single_scalers[0].initial_keys
    self._params = params
    self._experiments = experiments[0]
    self._Ih_table = JointIhTable(self.single_scalers)

  @abc.abstractmethod
  def calc_absorption_constraint(self, apm):
    'abstract method for combining the absorption constraints for multiple scalers'
    pass

  @staticmethod
  def calc_multi_absorption_constraint(apm, scalers):
    'method only called in aimless scaling'
    R = flex.double([])
    G = flex.double([])
    scaler_ids = [scaler.id_ for scaler in scalers]
    #if 'aimless' in scaler_ids:
    for i, scaler in enumerate(scalers):
      R.extend(scaler.calc_absorption_constraint(apm.apm_list[i])[0])
      G.extend(scaler.calc_absorption_constraint(apm.apm_list[i])[1])
    #else:
    #  for i, scaler in enumerate(scalers):

    return (R, G)

  @abc.abstractmethod
  def join_multiple_datasets(self):
    'abstract method for combining all datasets into a single reflectiont table'
    pass

  def join_datasets_from_scalers(self, scalers):
    '''method to create a joint reflection table from single scalers.
    Anticipated to be called from join_multiple_datasets'''
    joined_reflections = flex.reflection_table()
    joined_outliers = flex.reflection_table()
    for scaler in scalers:
      joined_reflections.extend(scaler.reflection_table)
      joined_outliers.extend(scaler.outlier_table)
    self._outlier_table = joined_outliers
    miller_set = miller.set(crystal.symmetry(
      space_group=scalers[0].experiments.crystal.get_space_group()),
      indices=joined_reflections['asu_miller_index'], anomalous_flag=False)
    permuted = miller_set.sort_permutation(by_value='packed_indices')
    self._reflection_table = joined_reflections.select(permuted)
    #does it make sense to do the following for targetscaler as well? do we want this?
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

class MultiScaler(MultiScalerBase):
  '''
  Scaler for multiple datasets - takes in params, experiments and
  a list of SingleScalers.
  '''

  id_ = 'multi'

  def __init__(self, params, experiments, single_scalers):
    logger.info('\nInitialising a MultiScaler instance. \n')
    super(MultiScaler, self).__init__(params, experiments, single_scalers)
    logger.info('Completed initialisation of MultiScaler. \n' + '*'*40 + '\n')

  def calc_absorption_constraint(self, apm):
    return super(MultiScaler, MultiScaler).calc_multi_absorption_constraint(
      apm, self.single_scalers)

  def update_error_model(self):
    for scaler in self.single_scalers:
      scaler.update_error_model()
    self._Ih_table = JointIhTable(self.single_scalers)

  def update_for_minimisation(self, apm):
    '''update the scale factors and Ih for the next iteration of minimisation,
    update the x values from the amp to the individual apms, as this is where
    basis functions, target functions etc get access to the parameters.'''
    for i, single_apm in enumerate(apm.apm_list):
      single_apm.x = apm.select_parameters(i)
    apm.derivatives = sparse.matrix(self.Ih_table.size, apm.n_active_params)
    for i, scaler in enumerate(self.single_scalers):
      basis_fn = scaler.get_basis_function(apm.apm_list[i])
      scaler.Ih_table.inverse_scale_factors = basis_fn[0]
      if basis_fn[1]:
        expanded = basis_fn[1].transpose() * self.Ih_table.h_index_expand_list[i]
        apm.derivatives.assign_block(expanded.transpose(), 0, apm.apm_data[i]['start_idx'])
    self.Ih_table.calc_Ih()

  def expand_scales_to_all_reflections(self, caller=None):
    for scaler in self.single_scalers:
      scaler.expand_scales_to_all_reflections(caller=self)

  def join_multiple_datasets(self):
    '''method to create a joint reflection table'''
    super(MultiScaler, self).join_datasets_from_scalers(self.single_scalers)

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

class TargetScaler(MultiScalerBase):
  '''
  Target Scaler for scaling one dataset against already scaled data - takes in
  params, lists of scaled and unscaled experiments, a list of already scaled
  SingleScalers and a list of unscaled reflections.
  '''

  id_ = 'target'

  def __init__(self, params, scaled_experiments, scaled_scalers,
    unscaled_experiments, unscaled_scalers):
    logger.info('\nInitialising a TargetScaler instance. \n')
    super(TargetScaler, self).__init__(params, scaled_experiments, scaled_scalers)
    self.unscaled_scalers = unscaled_scalers
    self._experiments = unscaled_experiments[0]
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
    logger.info('Completed initialisation of TargetScaler. \n' + '*'*40 + '\n')

  def get_target_function(self, apm):
    '''override the target function method for fixed Ih'''
    R = 0.0
    G = flex.double([])
    for i, scaler in enumerate(self.unscaled_scalers):
      (Ri, Gi) = target_function_fixedIh(scaler, apm.apm_list[i]).return_targets()
      R += Ri
      G.extend(Gi)
    return R, G

  def update_for_minimisation(self, apm):
    '''update the scale factors and Ih for the next iteration of minimisation'''
    for i, single_apm in enumerate(apm.apm_list):
      single_apm.x = apm.select_parameters(i)
    for i, scaler in enumerate(self.unscaled_scalers):
      basis_fn = scaler.get_basis_function(apm.apm_list[i])
      apm.apm_list[i].derivatives = basis_fn[1]
      scaler.Ih_table.inverse_scale_factors = basis_fn[0]

  def calc_absorption_constraint(self, apm):
    return super(TargetScaler, TargetScaler).calc_multi_absorption_constraint(
      apm, self.unscaled_scalers)

  def expand_scales_to_all_reflections(self, caller=None):
    for scaler in self.unscaled_scalers:
      scaler.expand_scales_to_all_reflections(caller=self)

  def calc_merging_statistics(self):
    results = []
    scaled_ids = []
    for scaler in self.unscaled_scalers:
      result, scaled_id = scaler.calc_merging_statistics()
      results.append(result[0])
      scaled_ids.append(scaled_id)
    return (results, scaled_ids)

  def join_multiple_datasets(self):
    '''method to create a joint reflection table'''
    scalers = []
    scalers.extend(self.single_scalers)
    scalers.extend(self.unscaled_scalers)
    super(TargetScaler, self).join_datasets_from_scalers(scalers)
