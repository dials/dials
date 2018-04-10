"""
This module defines a set of 'Scalers'. These act to initialise and connect
various parts of the scaling algorithm and datastructures such as the Ih_table,
basis_function etc, and present a united interface to the main scale.py
script. A SingleScaler is defined, for scaling of a single dataset, and a
MultiScaler is defined for scaling multiple datasets simultaneously.
A TargetScaler is used for targeted scaling.
"""

import abc
import logging
from dials.array_family import flex
from cctbx import miller, crystal
from scitbx import sparse
from dials_scaling_helpers_ext import row_multiply
import iotbx.merging_statistics
from libtbx.containers import OrderedSet
from libtbx.table_utils import simple_table
from dials.algorithms.scaling.basis_functions import basis_function
from dials.algorithms.scaling.scaling_utilities import (
  reject_outliers, calculate_wilson_outliers, calc_normE2)
#from dials.algorithms.scaling.restraints.scaling_restraints import ScalingRestraints
from dials.algorithms.scaling.reflection_weighting import Weighting
from dials.algorithms.scaling.Ih_table import SingleIhTable,\
  JointIhTable#, IhTableBase

from dials_scratch_scaling_ext import calc_sigmasq as cpp_calc_sigmasq
logger = logging.getLogger('dials')

class ScalingRestraints(object):
  """Scaling Restraints class."""
  def __init__(self, apm):
    self.apm = apm

  def compute_restraints_residuals_jacobian(self):
    """Calculate restraints for jacobian."""
    restraints = None
    restr = self.calculate_restraints()
    if restr:
      resid_restr = restr[0] # list
      surface_weights = self.apm.components['absorption']['object'].parameter_restraints
      n_abs_params = restr[0].size()
      n_tot_params = self.apm.n_active_params
      jacobian = sparse.matrix(n_abs_params, n_tot_params)
      offset = n_tot_params - n_abs_params
      for i in range(n_abs_params):
        jacobian[i, offset+i] = -1.0 * (surface_weights[i]**0.5)
      restraints = [resid_restr, jacobian, flex.double(n_abs_params, 1.0)]
    return restraints

  def calculate_restraints(self):
    """Calculate restraints for the scaling model."""
    residuals = flex.double([])
    gradient_vector = flex.double([])
    for comp in self.apm.components.itervalues():
      resid = comp['object'].calculate_restraints()
      if resid:
        gradient_vector.extend(resid[1])
        residuals.extend(resid[0])
      else:
        gradient_vector.extend(flex.double(comp['n_params'], 0.0))
    return [residuals, gradient_vector]

class ScalerBase(object):
  """Base class for all Scalers (single and multiple)."""

  __metaclass__ = abc.ABCMeta

  def __init__(self):
    self._experiments = None
    self._params = None
    self._reflection_table = []
    self._Ih_table = None
    self._initial_keys = []

  @property
  def Ih_table(self):
    """"A sorted reflection table, with additional methods for performing
    sums over equivalent reflections."""
    return self._Ih_table

  @Ih_table.setter
  def Ih_table(self, new_Ih_table):
    assert new_Ih_table.id_ == 'IhTableBase'
    self._Ih_table = new_Ih_table

  @property
  def experiments(self):
    """The experiment object associated with the dataset."""
    return self._experiments

  @property
  def reflection_table(self):
    """The reflection table of the datatset."""
    return self._reflection_table

  @property
  def params(self):
    """The params phil scope."""
    return self._params

  @property
  def initial_keys(self):
    """A list of initial reflection table keys."""
    return self._initial_keys

  def clean_reflection_table(self):
    """Remove additional added columns that are not required for output."""
    self._initial_keys.append('inverse_scale_factor')
    self._initial_keys.append('inverse_scale_factor_variance')
    self._initial_keys.append('Ih_values')
    for key in self.reflection_table.keys():
      if not key in self._initial_keys:
        del self._reflection_table[key]

  @abc.abstractmethod
  def update_for_minimisation(self, apm, curvatures=False):
    """Update the scale factors and Ih for the next minimisation iteration."""
    pass

  @abc.abstractmethod
  def expand_scales_to_all_reflections(self, caller=None):
    """Expand scales from a subset to all reflections."""
    pass

  def get_basis_function(self, apm, curvatures=False):
    """Call the basis function"""
    return basis_function(self, apm, curvatures).return_basis()

  @staticmethod
  def _map_indices_to_asu(reflection_table, experiments, params):
    """Map the miller index to the asu and use to sort the reflection table."""
    u_c = experiments.crystal.get_unit_cell().parameters()
    if params.scaling_options.space_group:
      sg_from_file = experiments.crystal.get_space_group().info()
      s_g_symbol = params.scaling_options.space_group
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
    permuted = (miller_set.map_to_asu()).sort_permutation(
      by_value='packed_indices')
    reflection_table = reflection_table.select(permuted)
    return reflection_table

  @classmethod
  def _scaling_subset(cls, reflection_table, params, error_model_params=None):
    """Select reflections with non-zero weight and update scale weights."""
    weights_for_scaling = cls._update_weights_for_scaling(reflection_table,
      params, error_model_params=error_model_params)
    sel = weights_for_scaling.weights > 0.0
    sel1 = reflection_table['Esq'] > params.reflection_selection.E2_min
    sel2 = reflection_table['Esq'] < params.reflection_selection.E2_max
    selection = sel & sel1 & sel2
    weights_for_scaling.weights.set_selected(~selection, 0.0)
    reflections_for_scaling = reflection_table.select(selection)
    #weights_for_scaling.weights = weights_for_scaling.weights.select(selection)
    msg = ('{0} reflections were selected for scale factor determination {sep}'
      'out of {5} reflections. This was based on selection criteria of {sep}'
      'E2_min = {1}, E2_max = {2}, Isigma_min = {3}, dmin = {4}. {sep}').format(
      reflections_for_scaling.size(), params.reflection_selection.E2_min,
      params.reflection_selection.E2_max, params.reflection_selection.Isigma_min,
      params.reflection_selection.d_min, reflection_table.size(), sep='\n')
    logger.info(msg)
    #return reflections_for_scaling, weights_for_scaling, selection
    return reflection_table, weights_for_scaling, selection

  @staticmethod
  def _update_weights_for_scaling(reflection_table, params,
    weights_filter=True, error_model_params=None):
    """Set the weights of each reflection to be used in scaling."""
    weights_for_scaling = Weighting(reflection_table)
    if weights_filter:
      weights_for_scaling.apply_Isigma_cutoff(reflection_table,
        params.reflection_selection.Isigma_min)
      weights_for_scaling.apply_dmin_cutoff(reflection_table,
        params.reflection_selection.d_min)
    if error_model_params:
      weights_for_scaling.apply_error_model(reflection_table,
        error_model_params)
    return weights_for_scaling

  def calc_merging_statistics(self):
    """Calculate the merging stats and return these with the dataset id."""
    u_c = self.experiments.crystal.get_unit_cell().parameters()
    bad_refl_sel = self.reflection_table.get_flags(
      self.reflection_table.flags.bad_for_scaling, all=False)
    reflections = self.reflection_table.select(~bad_refl_sel)
    if self.params.scaling_options.space_group:
      s_g_symbol = self.params.scaling_options.space_group
      crystal_symmetry = crystal.symmetry(unit_cell=u_c,
        space_group_symbol=s_g_symbol)
    else:
      s_g = self.experiments.crystal.get_space_group()
      crystal_symmetry = crystal.symmetry(unit_cell=u_c,
        space_group=s_g)
    miller_set = miller.set(crystal_symmetry=crystal_symmetry,
      indices=reflections['miller_index'], anomalous_flag=False)
    scaled_intensities = (reflections['intensity']/
      reflections['inverse_scale_factor'])
    sigmas = ((reflections['variance']**0.5)/
      reflections['inverse_scale_factor'])
    i_obs = miller.array(miller_set, data=scaled_intensities)
    i_obs.set_observation_type_xray_intensity()
    i_obs.set_sigmas(sigmas)
    scaled_ids = list(set(self.reflection_table['id']))
    result = iotbx.merging_statistics.dataset_statistics(
      i_obs=i_obs, n_bins=20, anomalous=False, sigma_filtering=None,
      use_internal_variance=True, eliminate_sys_absent=False)
      #eliminate_sys_absent=False, sigma_filtering=None)
    return ([result], scaled_ids[0])

class SingleScalerBase(ScalerBase):
  """
  Parent class for single-dataset Scalers, containing a standard
  setup routine for the reflection_table - takes in params, experiment
  and reflection.
  """

  id_ = 'single'

  def __init__(self, params, experiment, reflection, scaled_id=0):
    logger.info('Configuring a Scaler for a single dataset. \n')
    logger.info('The dataset id assigned to this reflection table is %s, \n'
      'the scaling model type being applied is %s. \n', scaled_id,
      experiment.scaling_model.id_)
    super(SingleScalerBase, self).__init__()
    self._experiments = experiment
    self._params = params
    n_model_params = sum([val.n_params for val in self.components.itervalues()])
    self._var_cov = sparse.matrix(n_model_params, n_model_params)
    # Rename the id for now so that we have unique dataset ids
    reflection['id'] = flex.int([scaled_id] * reflection.size())
    self._initial_keys = [key for key in reflection.keys()]
    # Choose intensities, map to asu, assign unique refl. index
    reflection_table = self._reflection_table_setup(self._initial_keys,
      reflection) #flag not integrated, bad d, bad partiality
    reflection_table = self._select_optimal_intensities(reflection_table,
      self.params) #flag bad variance,
    reflection_table = self._map_indices_to_asu(reflection_table,
      self.experiments, self.params)
    # Calculate values for later filtering, but don't filter here!!!
    reflection_table = calc_normE2(reflection_table, self.experiments)
    if 'centric_flag' in reflection_table: #If E2 values have been calculated
      reflection_table = calculate_wilson_outliers(reflection_table)
    if self.params.scaling_options.reject_outliers:
      reflection_table = self.round_of_outlier_rejection(reflection_table)
    self._reflection_table = reflection_table
    self._configure_reflection_table()
    self.select_reflections_for_scaling()
    logger.info('Completed configuration of Scaler. \n\n' + '='*80 + '\n')

  @property
  def components(self):
    """Shortcut to scaling model components."""
    return self.experiments.scaling_model.components

  @property
  def consecutive_refinement_order(self):
    """Shortcut to scaling model scaling order."""
    return self.experiments.scaling_model.consecutive_refinement_order

  @property
  def var_cov_matrix(self):
    return self._var_cov

  def update_var_cov(self, apm):
    """
    Update the full parameter variance covariance matrix after a refinement.

    If all parameters have been refined, then the full var_cov matrix can be set.
    Else one must select subblocks for pairs of parameters and assign these into
    the full var_cov matrix, taking care to out these in the correct position.
    This is applicable if only some parameters have been refined in this cycle.
    """
    var_cov_list = apm.var_cov_matrix #values are passed as a list from refinery
    if int(var_cov_list.size()**0.5) == self.var_cov_matrix.n_rows:
      self._var_cov.assign_block(var_cov_list.matrix_copy_block(0, 0,
        apm.n_active_params, apm.n_active_params), 0, 0)
    else: #need to set part of the var_cov matrix e.g. if only refined some params
      #first work out the order in self._var_cov
      cumul_pos_dict = {}
      n_cumul_params = 0
      for name, component in self.components.iteritems():
        cumul_pos_dict[name] = n_cumul_params
        n_cumul_params += component.n_params
      #now get a var_cov_matrix subblock for pairs of parameters
      for name in apm.components_list:
        for name2 in apm.components_list:
          n_rows = apm.components[name]['n_params']
          n_cols = apm.components[name2]['n_params']
          start_row = apm.components[name]['start_idx']
          start_col = apm.components[name2]['start_idx']
          sub = var_cov_list.matrix_copy_block(start_row, start_col, n_rows,
            n_cols)
          #now set this block into correct location in overall var_cov
          self._var_cov.assign_block(sub, cumul_pos_dict[name],
            cumul_pos_dict[name2])

  def update_for_minimisation(self, apm, curvatures=False):
    """Update the scale factors and Ih for the next minimisation iteration."""
    basis_fn = self.get_basis_function(apm, curvatures=curvatures)
    apm.derivatives = basis_fn[1]
    if curvatures:
      apm.curvatures = basis_fn[2]
    self.Ih_table.inverse_scale_factors = basis_fn[0]
    self.Ih_table.calc_Ih()

  @staticmethod
  def _reflection_table_setup(initial_keys, reflections):
    """Initial filter to select integrated reflections."""
    mask = ~reflections.get_flags(reflections.flags.integrated)
    d_mask = reflections['d'] <= 0.0
    partials_mask = reflections['partiality'] < 0.95
    reflections.set_flags(mask or partials_mask or d_mask,
      reflections.flags.excluded_for_scaling)
    print('%s reflections not suitable for scaling (low partiality,\n'
      'not integrated etc).\n' % reflections.get_flags(
      reflections.flags.excluded_for_scaling).count(True))
    if not 'inverse_scale_factor' in initial_keys:
      reflections['inverse_scale_factor'] = (flex.double(reflections.size(), 1.0))
    return reflections

  @classmethod
  def _select_optimal_intensities(cls, reflection_table, params):
    """Choose which intensities to use for scaling."""
    if (params.scaling_options.integration_method == 'sum' or
        params.scaling_options.integration_method == 'prf'):
      intstr = params.scaling_options.integration_method
      if 'dqe' in reflection_table:
        conversion = reflection_table['lp'] / reflection_table['dqe']
      else:
        conversion = reflection_table['lp']
      reflection_table['intensity'] = (
        reflection_table['intensity.'+intstr+'.value'] * conversion)
      reflection_table['variance'] = (
        reflection_table['intensity.'+intstr+'.variance'] * (conversion**2))
      logger.info('%s intensity values will be used for scaling. \n',
        'Profile fitted' if intstr == 'prf' else 'Summation integrated')
    #perform a combined prf/sum in a similar fashion to aimless
    elif params.scaling_options.integration_method == 'combine':
      int_prf = reflection_table['intensity.prf.value'] * conversion
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
      logger.info(msg)
    variance_mask = reflection_table['variance'] <= 0.0
    reflection_table.set_flags(variance_mask,
      reflection_table.flags.excluded_for_scaling)
    return reflection_table

  def expand_scales_to_all_reflections(self, caller=None, calc_cov=True):
    self._reflection_table['inverse_scale_factor'] = flex.double(
      self.reflection_table.size(), 1.0)
    self._reflection_table['inverse_scale_factor_variance'] = flex.double(
      self.reflection_table.size(), 0.0)
    sel = self.reflection_table.get_flags(
      self.reflection_table.flags.bad_for_scaling, all=False)
    scaled_sel = ~sel
    scaled_isel = scaled_sel.iselection()
    scaled_reflections = self._reflection_table.select(scaled_sel)
    scaled_reflections['inverse_scale_factor'] = flex.double(
      scaled_reflections.size(), 1.0)
    for component in self.components.itervalues():
      component.update_reflection_data(self.reflection_table, scaled_sel)
      component.calculate_scales_and_derivatives()
      scaled_reflections['inverse_scale_factor'] *= component.inverse_scales
    logger.info('Scale factors determined during minimisation have now been\n'
      'applied to all reflections for dataset %s.\n',
      self.reflection_table['id'][0])
    if self.var_cov_matrix and calc_cov:
      scaled_reflections['inverse_scale_factor_variance'] = calc_sf_variances(
        self.components, self._var_cov)
    self.reflection_table['inverse_scale_factor'].set_selected(scaled_isel,
      scaled_reflections['inverse_scale_factor'])
    self.reflection_table['inverse_scale_factor_variance'].set_selected(
      scaled_isel, scaled_reflections['inverse_scale_factor_variance'])
    if (self.params.scaling_options.reject_outliers and
      not isinstance(caller, MultiScalerBase)):
      self._reflection_table = self.round_of_outlier_rejection(
        self._reflection_table)
    if self.params.weighting.optimise_error_model:
      self.Ih_table = SingleIhTable(self._reflection_table)

  def update_error_model(self, error_model_params):
    """Apply a correction to try to improve the error estimate."""
    self.Ih_table.update_error_model(error_model_params)
    self.experiments.scaling_model.set_error_model(list(error_model_params))

  def compute_restraints_residuals_jacobian(self, apm):
    """Calculate a restraint for the jacobian."""
    restraints = None
    if self.experiments.scaling_model.id_ == 'physical':
      restraints = ScalingRestraints(apm).compute_restraints_residuals_jacobian()
    return restraints

  def calculate_restraints(self, apm):
    """Calculate a restraint for the residual and gradient."""
    restraints = None
    if self.experiments.scaling_model.id_ == 'physical':
      restraints = ScalingRestraints(apm).calculate_restraints()
    return restraints

  def round_of_outlier_rejection(self, reflection_table):
    """calculate outliers from the reflections in the Ih_table,
    and use these to filter the reflection table and Ih_table."""
    reflection_table = reject_outliers(reflection_table,
      self.params.scaling_options.outlier_zmax)
    msg = ('A round of outlier rejection has been performed, in total {0} {sep}'
        'outliers have now been identified. {sep}'.format(
        reflection_table.get_flags(reflection_table.flags.outlier_in_scaling
        ).count(True), sep='\n'))
    logger.info(msg)
    return reflection_table

  def apply_error_model_to_variances(self):
    error_model = self.experiments.scaling_model.configdict['error_model_parameters']
    self.reflection_table['variance'] = error_model[0] * (self.reflection_table['variance']
      + ((error_model[1] * self.reflection_table['intensity'])**2))**0.5
    logger.info('The error model determined has been applied to the variances')

  def apply_selection_to_SFs(self, sel):
    """Updates data from within current SF objects using the given selection.
       Required for targeted scaling."""
    for component in self.components.itervalues():
      component.update_reflection_data(self.reflection_table, sel)

  def select_reflections_for_scaling(self):
    """Select a subset of reflections, create and Ih table and update the
    model components."""
    (refl_for_scaling, w_for_scaling, selection) = (
      self._scaling_subset(self.reflection_table, self.params))
    self._Ih_table = SingleIhTable(refl_for_scaling, w_for_scaling.weights)
    #if self.params.weighting.tukey_biweighting: #not working for now, FIX
    #  self.Ih_table.apply_tukey_biweighting()
    for component in self.components.itervalues():
      component.update_reflection_data(self.reflection_table, selection)

  def _configure_reflection_table(self):
    """Calculate requried quantities"""
    self._reflection_table = self.experiments.scaling_model.configure_reflection_table(
      self._reflection_table, self.experiments, self.params)
    rows = []
    for key, val in self.components.iteritems():
      rows.append([key, str(val.n_params)])
    st = simple_table(rows, ['correction', 'n_parameters'])
    logger.info('The following corrections will be applied to this dataset: \n')
    logger.info(st.format())

class MultiScalerBase(ScalerBase):
  '''Base class for Scalers handling multiple datasets'''
  def __init__(self, params, experiments, single_scalers):
    '''initialise from a list of single scalers'''
    super(MultiScalerBase, self).__init__()
    self.single_scalers = single_scalers
    self._initial_keys = self.single_scalers[0].initial_keys
    self._params = params
    self._experiments = experiments[0] #what should this be set to - where exactly is used?
    logger.info('Determining symmetry equivalent reflections across datasets.\n')
    self._Ih_table = JointIhTable(self.single_scalers)

  @abc.abstractmethod
  def calculate_restraints(self, apm):
    'abstract method for combining the absorption restraints for multiple scalers'
    pass

  @staticmethod
  def calc_multi_absorption_restraint(apm, scalers):
    'method only called in physical scaling'
    restraints = None
    scaler_ids = [scaler.experiments.scaling_model.id_ for scaler in scalers]
    if 'physical' in scaler_ids:
      R = flex.double([])
      G = flex.double([])
      for i, scaler in enumerate(scalers):
        restr = scaler.calculate_restraints(apm.apm_list[i])
        if restr:
          R.extend(restr[0])
          G.extend(restr[1])
        else:
          G.extend(flex.double(apm.apm_list[i].n_active_params, 0.0))
      restraints = [R, G]
    return restraints

  @abc.abstractmethod
  def compute_restraints_residuals_jacobian(self, apm):
    'abstract method for combining the absorption restraints for multiple scalers'
    pass

  @staticmethod
  def compute_multi_restraints_residuals_jacobian(apm, scalers):
    '''method to calculate absorption restraint for multiple scalers'''
    restraints = None
    scaler_ids = [scaler.experiments.scaling_model.id_ for scaler in scalers]
    if 'physical' in scaler_ids:
      R = flex.double([])
      jacobian_list = []
      scaler_list = []
      for i, scaler in enumerate(scalers):
        restr = scaler.compute_restraints_residuals_jacobian(apm.apm_list[i])
        if restr:
          R.extend(restr[0])
          jacobian_list.append(restr[1])
          scaler_list.append(i)
      n_restraints = R.size()
      jacobian = sparse.matrix(n_restraints, apm.n_active_params)
      n_row_cumul = 0
      for j, jac in enumerate(jacobian_list):
        scaler_no = scaler_list[j]
        jacobian.assign_block(jac, n_row_cumul, apm.apm_data[scaler_no]['start_idx'])
        n_row_cumul += jac.n_rows
      restraints = [R, jacobian, flex.double(R.size(), 1.0)]
    return restraints

  @abc.abstractmethod
  def join_multiple_datasets(self):
    'abstract method for combining all datasets into a single reflectiont table'
    pass

  def join_datasets_from_scalers(self, scalers):
    '''method to create a joint reflection table from single scalers.
    Anticipated to be called from join_multiple_datasets'''
    joined_reflections = flex.reflection_table()
    for scaler in scalers:
      joined_reflections.extend(scaler.reflection_table)
    miller_set = miller.set(crystal.symmetry(
      space_group=scalers[0].experiments.crystal.get_space_group()),
      indices=joined_reflections['asu_miller_index'], anomalous_flag=False)
    permuted = miller_set.sort_permutation(by_value='packed_indices')
    self._reflection_table = joined_reflections.select(permuted)
    if self.params.scaling_options.reject_outliers:
      self._reflection_table = self.round_of_outlier_rejection(
        self._reflection_table)

  def round_of_outlier_rejection(self, reflection_table):
    """calculate outliers from the reflections in the Ih_table,
    and use these to filter the reflection table and Ih_table."""
    reflection_table = reject_outliers(reflection_table,
      self.params.scaling_options.outlier_zmax)
    msg = ('Combined outlier rejection has been performed across all datasets, {sep}'
      'in total {0} outliers have now been identified. {sep}'.format(
        reflection_table.get_flags(reflection_table.flags.outlier_in_scaling
        ).count(True), sep='\n'))
    logger.info(msg)
    return reflection_table

  def apply_error_model_to_variances(self):
    """Update variances of individual reflection tables."""
    for scaler in self.single_scalers:
      scaler.apply_error_model_to_variances()

  def select_reflections_for_scaling(self):
    for scaler in self.single_scalers:
      scaler.select_reflections_for_scaling()

class MultiScaler(MultiScalerBase):
  """
  Scaler for multiple datasets - takes in params, experiments and
  a list of SingleScalers.
  """

  id_ = 'multi'

  def __init__(self, params, experiments, single_scalers):
    logger.info('Configuring a MultiScaler to handle the individual Scalers. \n')
    super(MultiScaler, self).__init__(params, experiments, single_scalers)
    logger.info('Completed configuration of MultiScaler. \n\n' + '='*80 + '\n')

  def calculate_restraints(self, apm):
    return super(MultiScaler, MultiScaler).calc_multi_absorption_restraint(
      apm, self.single_scalers)

  def compute_restraints_residuals_jacobian(self, apm):
    return super(MultiScaler, MultiScaler).compute_multi_restraints_residuals_jacobian(
      apm, self.single_scalers)

  def update_error_model(self, error_model_params):
    self.Ih_table.update_error_model(error_model_params)

  def update_for_minimisation(self, apm, curvatures=False):
    '''update the scale factors and Ih for the next iteration of minimisation,
    update the x values from the amp to the individual apms, as this is where
    basis functions, target functions etc get access to the parameters.'''
    for i, single_apm in enumerate(apm.apm_list):
      single_apm.x = apm.select_parameters(i)
    apm.derivatives = sparse.matrix(self.Ih_table.size, apm.n_active_params)
    start_row_no = 0
    for i, scaler in enumerate(self.single_scalers):
      basis_fn = scaler.get_basis_function(apm.apm_list[i])
      scaler.Ih_table.inverse_scale_factors = basis_fn[0]
      if basis_fn[1]:
        apm.derivatives.assign_block(basis_fn[1], start_row_no,
          apm.apm_data[i]['start_idx'])
        start_row_no += basis_fn[1].n_rows
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
  """
  Target Scaler for scaling one dataset against already scaled data - takes in
  params, lists of scaled and unscaled experiments, a list of already scaled
  SingleScalers and a list of unscaled reflections.
  """

  id_ = 'target'

  def __init__(self, params, scaled_experiments, scaled_scalers,
    unscaled_experiments, unscaled_scalers):
    logger.info('\nInitialising a TargetScaler instance. \n')
    super(TargetScaler, self).__init__(params, scaled_experiments, scaled_scalers)
    self.unscaled_scalers = unscaled_scalers
    self._initial_keys = self.unscaled_scalers[0].initial_keys #needed for
    # scaling against calculated Is
    self._experiments = unscaled_experiments[0]
    target_Ih_table = self.Ih_table
    target_asu_Ih_dict = dict(zip(target_Ih_table.asu_miller_index,
      target_Ih_table.Ih_values))
    for scaler in unscaled_scalers:
      scaler.Ih_table.Ih_table['Ih_values'] = flex.double(
        scaler.Ih_table.size, 0.0) # set to zero to allow selection below
      location_in_unscaled_array = 0
      for j, miller_idx in enumerate(OrderedSet(scaler.Ih_table.asu_miller_index)):
        n_in_group = scaler.Ih_table.h_index_matrix.col(j).non_zeroes
        if miller_idx in target_asu_Ih_dict:
          i = location_in_unscaled_array
          Ih = flex.double(n_in_group, target_asu_Ih_dict[miller_idx])
          scaler.Ih_table.Ih_values.set_selected(
            flex.size_t(range(i, i + n_in_group)), Ih)
        location_in_unscaled_array += n_in_group
      sel = scaler.Ih_table.Ih_values != 0.0
      scaler.Ih_table = scaler.Ih_table.select(sel)
      scaler.apply_selection_to_SFs(scaler.Ih_table.nonzero_weights)
    logger.info('Completed initialisation of TargetScaler. \n' + '*'*40 + '\n')

  def update_for_minimisation(self, apm, curvatures=False):
    """Update the scale factors and Ih for the next iteration of minimisation."""
    for i, single_apm in enumerate(apm.apm_list):
      single_apm.x = apm.select_parameters(i)
    for i, scaler in enumerate(self.unscaled_scalers):
      basis_fn = scaler.get_basis_function(apm.apm_list[i])
      apm.apm_list[i].derivatives = basis_fn[1]
      scaler.Ih_table.inverse_scale_factors = basis_fn[0]

  def calculate_restraints(self, apm):
    return super(TargetScaler, TargetScaler).calc_multi_absorption_restraint(
      apm, self.unscaled_scalers)

  def compute_restraints_residuals_jacobian(self, apm):
    return super(TargetScaler, TargetScaler).compute_multi_restraints_residuals_jacobian(
      apm, self.unscaled_scalers)

  def expand_scales_to_all_reflections(self, caller=None, calc_cov=True):
    for scaler in self.unscaled_scalers:
      scaler.expand_scales_to_all_reflections(caller=self, calc_cov=calc_cov)

  def calc_merging_statistics(self):
    results = []
    scaled_ids = []
    for scaler in self.unscaled_scalers:
      result, scaled_id = scaler.calc_merging_statistics()
      results.append(result[0])
      scaled_ids.append(scaled_id)
    joint_result, _ = super(TargetScaler, self).calc_merging_statistics()
    results.append(joint_result[0])
    scaled_ids.append('x')
    return (results, scaled_ids)

  def join_multiple_datasets(self):
    '''method to create a joint reflection table'''
    scalers = []
    if not self.params.scaling_options.target_intensities:
      scalers.extend(self.single_scalers)
    scalers.extend(self.unscaled_scalers)
    super(TargetScaler, self).join_datasets_from_scalers(scalers)

  def select_reflections_for_scaling(self):
    for scaler in self.unscaled_scalers:
      scaler.select_reflections_for_scaling()


class NullScaler(ScalerBase):
  """
  Null scaler to allow targeted scaling against calculated intensities.
  """

  id_ = 'null'

  def __init__(self, params, experiment, reflection, scaled_id=0):
    super(NullScaler, self).__init__()
    self._experiments = experiment
    self._params = params
    self._scaled_id = scaled_id
    self._reflection_table = self._map_indices_to_asu(reflection, experiment,
      params)
    self._initial_keys = [key for key in self._reflection_table.keys()]
    self._reflection_table['intensity'] = self._reflection_table[
      'intensity.calculated.value']
    n_refl = self._reflection_table.size()
    self._reflection_table['inverse_scale_factor'] = flex.double(n_refl, 1.0)
    self._reflection_table['variance'] = flex.double(n_refl, 1.0)
    weights = self._reflection_table['intensity.calculated.value']
    self._reflection_table.set_flags(flex.bool([False]*n_refl),
      self._reflection_table.flags.excluded_for_scaling)
    self._Ih_table = SingleIhTable(self._reflection_table, weights)
    logger.info('NullScaler contains %s reflections', n_refl)
    logger.info('Completed configuration of NullScaler. \n\n' + '='*80 + '\n')

  def expand_scales_to_all_reflections(self, caller=None):
    pass

  def update_for_minimisation(self, apm, curvatures=False):
    '''update the scale factors and Ih for the next iteration of minimisation'''
    pass

def calc_sf_variances(components, var_cov):
  """Use the parameter var_cov matrix to calculate the variances of the
  inverse scales."""
  n_param = 0
  for component in components:
    n_param += components[component].n_params
    n_refl = components[component].inverse_scales.size() #should all be same
  jacobian = sparse.matrix(n_refl, n_param)
  n_cumulative_param = 0
  for component in components:
    block = components[component].derivatives
    n_param = components[component].n_params
    for component_2 in components:
      if component_2 != component:
        block = row_multiply(block, components[component].inverse_scales)
    jacobian.assign_block(block, 0, n_cumulative_param)
    n_cumulative_param += n_param
  logger.info('Calculating error estimates of inverse scale factors. \n')
  return cpp_calc_sigmasq(jacobian.transpose(), var_cov)
