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
from cctbx import crystal
from scitbx import sparse
from dials_scaling_helpers_ext import row_multiply
from libtbx.table_utils import simple_table
from dials.algorithms.scaling.basis_functions import basis_function
from dials.algorithms.scaling.scaling_utilities import (
  set_wilson_outliers, calc_normE2)
from dials.algorithms.scaling.outlier_rejection import reject_outliers
from dials.algorithms.scaling.Ih_table import SingleIhTable,\
  JointIhTable
from dials.algorithms.scaling.scaling_restraints import ScalingRestraints,\
  MultiScalingRestraints
from dials.algorithms.scaling.target_function import ScalingTarget,\
  ScalingTargetFixedIH
from dials.algorithms.scaling.scaling_refiner import scaling_refinery,\
  error_model_refinery
from dials.algorithms.scaling.error_model.error_model import \
  BasicErrorModel
from dials.algorithms.scaling.error_model.error_model_target import \
  ErrorModelTarget
from dials.algorithms.scaling.parameter_handler import create_apm_factory
from dials_scratch_scaling_ext import calc_sigmasq as cpp_calc_sigmasq
logger = logging.getLogger('dials')


class ScalerBase(object):
  """Base class for all Scalers (single and multiple)."""

  __metaclass__ = abc.ABCMeta

  def __init__(self):
    self._experiments = None
    self.space_group = None
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
    msg = ('Attempting to set a new Ih_table with an object not recognised as \n'
      'a valid Ih_table')
    assert hasattr(new_Ih_table, 'id_'), msg
    assert new_Ih_table.id_ == 'IhTableBase', msg
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
  def expand_scales_to_all_reflections(self, caller=None, calc_cov=False):
    """Expand scales from a subset to all reflections."""
    pass

  def _set_space_group(self):
    if self.params.scaling_options.space_group:
      sg_from_file = self.experiments.crystal.get_space_group().info()
      s_g_symbol = self.params.scaling_options.space_group
      crystal_symmetry = crystal.symmetry(space_group_symbol=s_g_symbol)
      self.space_group = crystal_symmetry.space_group()
      self.experiments.crystal.set_space_group(self.space_group)
      msg = ('WARNING: Manually overriding space group from {0} to {1}. {sep}'
        'If the reflection indexing in these space groups is different, {sep}'
        'bad things may happen!!! {sep}').format(sg_from_file, s_g_symbol, sep='\n')
      logger.info(msg)
    else:
      self.space_group = self.experiments.crystal.get_space_group()

  @classmethod
  def _scaling_subset(cls, reflection_table, params):
    """Select reflections with non-zero weight and update scale weights."""
    sel = ~reflection_table.get_flags(reflection_table.flags.bad_for_scaling,
      all=False)
    sel1 = reflection_table['Esq'] > params.reflection_selection.E2_min
    sel2 = reflection_table['Esq'] < params.reflection_selection.E2_max
    sel3 = reflection_table['intensity']/(reflection_table['variance']**0.5) > (
      params.reflection_selection.Isigma_min)
    sel4 = reflection_table['d'] > params.reflection_selection.d_min
    selection = sel & sel1 & sel2 & sel3 & sel4
    msg = ('{0} reflections were selected for scale factor determination {sep}'
      'out of {5} reflections. This was based on selection criteria of {sep}'
      'E2_min = {1}, E2_max = {2}, Isigma_min = {3}, dmin = {4}. {sep}').format(
      selection.count(True), params.reflection_selection.E2_min,
      params.reflection_selection.E2_max, params.reflection_selection.Isigma_min,
      params.reflection_selection.d_min, reflection_table.size(), sep='\n')
    logger.info(msg)
    return selection

  def perform_scaling(self, target_type=ScalingTarget, engine=None,
      max_iterations=None):
    """Minimise the scaling model"""
    apm_factory = create_apm_factory(self)
    for _ in range(apm_factory.n_cycles):
      apm = apm_factory.make_next_apm()
      if not engine:
        engine = self.params.scaling_refinery.engine
      if not max_iterations:
        max_iterations = self.params.scaling_refinery.max_iterations
      refinery = scaling_refinery(engine=engine, target=target_type(self, apm),
        prediction_parameterisation=apm, max_iterations=max_iterations)
      refinery.run()
      self = refinery.return_scaler()
      logger.info('\n'+'='*80+'\n')

  @abc.abstractmethod
  def perform_error_optimisation(self):
    """Optimise an error model."""
    pass


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
    self._set_space_group()
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
    # Calculate values for later filtering, but don't filter here!!!
    reflection_table = calc_normE2(reflection_table, self.experiments)
    #if 'centric_flag' in reflection_table: #If E2 values have been calculated
    #  reflection_table = set_wilson_outliers(reflection_table)
    if self.params.scaling_options.outlier_rejection != '0':
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
    """Link to consecutive refinement order for parameter manager."""
    return self.experiments.scaling_model.consecutive_refinement_order

  @property
  def var_cov_matrix(self):
    """The variance covariance matrix for the parameters."""
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
    basis_fn = basis_function(apm, curvatures).return_basis()
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
    partials_mask = reflections['partiality'] < 0.6
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
      if 'dqe' in reflection_table:
        conversion /= reflection_table['dqe']
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

  def expand_scales_to_all_reflections(self, caller=None, calc_cov=False):
    self._reflection_table['inverse_scale_factor'] = flex.double(
      self.reflection_table.size(), 1.0)
    self._reflection_table['inverse_scale_factor_variance'] = flex.double(
      self.reflection_table.size(), 0.0)
    scaled_sel = ~self.reflection_table.get_flags(
      self.reflection_table.flags.bad_for_scaling, all=False)
    scaled_isel = scaled_sel.iselection()
    scaled_invsf = self._reflection_table['inverse_scale_factor'].select(scaled_sel)
    for component in self.components.itervalues():
      component.update_reflection_data(self.reflection_table, scaled_sel)
      component.calculate_scales_and_derivatives()
      scaled_invsf *= component.inverse_scales
    self.reflection_table['inverse_scale_factor'].set_selected(scaled_isel,
      scaled_invsf)
    logger.info('Scale factors determined during minimisation have now been\n'
      'applied to all reflections for dataset %s.\n',
      self.reflection_table['id'][0])
    if self.var_cov_matrix and calc_cov:
      scaled_invsfvar = calc_sf_variances(self.components, self._var_cov)
      self.reflection_table['inverse_scale_factor_variance'].set_selected(
        scaled_isel, scaled_invsfvar)
    if (self.params.scaling_options.outlier_rejection != '0' and
      not isinstance(caller, MultiScalerBase)):
      self._reflection_table = self.round_of_outlier_rejection(
        self._reflection_table)
    #if self.params.weighting.optimise_error_model:
    #  self.Ih_table = SingleIhTable(self._reflection_table, self.space_group)

  def update_error_model(self, error_model_params):
    """Apply a correction to try to improve the error estimate."""
    self.Ih_table.update_error_model(error_model_params)
    self.experiments.scaling_model.set_error_model(list(error_model_params))

  def compute_restraints_residuals_jacobian(self, apm):
    """Calculate a restraint for the jacobian."""
    restraints = None
    if self.experiments.scaling_model.id_ == 'physical':
      restraints = ScalingRestraints(apm).calculate_jacobian_restraints()
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
    reflection_table = reject_outliers(reflection_table, self.space_group,
      self.params)
    msg = ('A round of outlier rejection has been performed, in total {0} {sep}'
        'outliers have now been identified. {sep}'.format(
        reflection_table.get_flags(reflection_table.flags.outlier_in_scaling
        ).count(True), sep='\n'))
    logger.info(msg)
    #exit()
    return reflection_table

  def apply_error_model_to_variances(self):
    """Apply an aimless-like error model to the variances."""
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
    selection = self._scaling_subset(self.reflection_table, self.params)
    self._Ih_table = SingleIhTable(self.reflection_table, self.space_group,
      selection, self.params.weighting.weighting_scheme)
    if self.params.weighting.error_model_params:
      self.update_error_model(self.params.weighting.error_model_params)
    for component in self.components.itervalues():
      component.update_reflection_data(self.reflection_table, selection)

  def perform_error_optimisation(self):
    self.Ih_table = SingleIhTable(self._reflection_table, self.space_group)
    refinery = error_model_refinery(engine='SimpleLBFGS',
      target=ErrorModelTarget(BasicErrorModel(self.Ih_table)),
      max_iterations=100)
    refinery.run()
    error_model = refinery.return_error_manager()
    self.update_error_model(error_model.refined_parameters)

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
    self.space_group = single_scalers[0].space_group
    self.active_scalers = None
    self._initial_keys = self.single_scalers[0].initial_keys
    self._params = params
    self._experiments = experiments[0]
    logger.info('Determining symmetry equivalent reflections across datasets.\n')
    self._Ih_table = JointIhTable([x.Ih_table for x in self.single_scalers])

  def calculate_restraints(self, apm):
    """Calculate a restraints residuals/gradient vector for multiple datasets."""
    scaler_ids = [scaler.experiments.scaling_model.id_ for scaler in
      self.active_scalers]
    if 'physical' in scaler_ids:
      #Ok for now, but make more general? But will be slow for many datasets?
      return MultiScalingRestraints(apm).calculate_restraints()
    return None

  def compute_restraints_residuals_jacobian(self, apm):
    """Calculate a restraints residuals/jacobian for multiple datasets."""
    scaler_ids = [scaler.experiments.scaling_model.id_ for scaler in
      self.active_scalers]
    if 'physical' in scaler_ids:
      return MultiScalingRestraints(apm).calculate_jacobian_restraints()
    return None

  def expand_scales_to_all_reflections(self, caller=None, calc_cov=False):
    for scaler in self.active_scalers:
      scaler.expand_scales_to_all_reflections(caller=self, calc_cov=calc_cov)

  @abc.abstractmethod
  def join_multiple_datasets(self):
    """Combine all datasets into a single reflection table."""
    pass

  def join_datasets_from_scalers(self, scalers):
    """Create a joint reflection table from single scalers.
    Anticipated to be called from join_multiple_datasets."""
    self._reflection_table = flex.reflection_table()
    for scaler in scalers:
      self._reflection_table.extend(scaler.reflection_table)
    if self.params.scaling_options.outlier_rejection != '0':
      self._reflection_table = self.round_of_outlier_rejection(
        self._reflection_table)

  def round_of_outlier_rejection(self, reflection_table):
    """calculate outliers from the reflections in the Ih_table,
    and use these to filter the reflection table and Ih_table."""
    reflection_table = reject_outliers(reflection_table, self.space_group,
      self.params)
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
    """Select reflections for scaling in individual scalers."""
    for scaler in self.active_scalers:
      scaler.select_reflections_for_scaling()

  def perform_error_optimisation(self):
    """Perform an optimisation of the sigma values."""
    for s_scaler in self.active_scalers:
      s_scaler.Ih_table = SingleIhTable(s_scaler.reflection_table,
        s_scaler.space_group)
      refinery = error_model_refinery(engine='SimpleLBFGS',
        target=ErrorModelTarget(BasicErrorModel(s_scaler.Ih_table)),
        max_iterations=100)
      refinery.run()
      error_model = refinery.return_error_manager()
      s_scaler.update_error_model(error_model.refined_parameters)
    self.Ih_table.update_weights_from_error_models()

class MultiScaler(MultiScalerBase):
  """
  Scaler for multiple datasets - takes in params, experiments and
  a list of SingleScalers.
  """

  id_ = 'multi'

  def __init__(self, params, experiments, single_scalers):
    logger.info('Configuring a MultiScaler to handle the individual Scalers. \n')
    super(MultiScaler, self).__init__(params, experiments, single_scalers)
    self.active_scalers = self.single_scalers
    logger.info('Completed configuration of MultiScaler. \n\n' + '='*80 + '\n')

  def update_error_model(self, error_model_params):
    """Update the error model in Ih table."""
    self.Ih_table.update_error_model(error_model_params)

  def update_for_minimisation(self, apm, curvatures=False):
    '''update the scale factors and Ih for the next iteration of minimisation,
    update the x values from the amp to the individual apms, as this is where
    basis functions, target functions etc get access to the parameters.'''
    apm.derivatives = sparse.matrix(self.Ih_table.size, apm.n_active_params)
    start_row_no = 0
    for i, scaler in enumerate(self.single_scalers):
      basis_fn = basis_function(apm.apm_list[i], curvatures).return_basis()
      scaler.Ih_table.inverse_scale_factors = basis_fn[0]
      if basis_fn[1]:
        apm.derivatives.assign_block(basis_fn[1], start_row_no,
          apm.apm_data[i]['start_idx'])
        start_row_no += basis_fn[1].n_rows
    self.Ih_table.calc_Ih()

  def join_multiple_datasets(self):
    '''method to create a joint reflection table'''
    super(MultiScaler, self).join_datasets_from_scalers(self.single_scalers)

class TargetScaler(MultiScalerBase):
  """
  Target Scaler for scaling one dataset against already scaled data - takes in
  params, lists of scaled and unscaled experiments, a list of already scaled
  SingleScalers and a list of unscaled reflections.
  """

  id_ = 'target'

  def __init__(self, params, scaled_experiments, scaled_scalers,
    unscaled_scalers):
    logger.info('\nInitialising a TargetScaler instance. \n')
    super(TargetScaler, self).__init__(params, scaled_experiments, scaled_scalers)
    self.unscaled_scalers = unscaled_scalers
    self.active_scalers = self.unscaled_scalers
    self._initial_keys = self.unscaled_scalers[0].initial_keys #needed for
    # scaling against calculated Is
    self.set_Ih_values_to_target()
    logger.info('Completed initialisation of TargetScaler. \n' + '*'*40 + '\n')

  def set_Ih_values_to_target(self):
    """Match equivalent reflections between individual unscaled datasets and
    set target Ixh values in the Ih_table for each dataset."""
    target_Ih_table = self.Ih_table
    for scaler in self.unscaled_scalers:
      scaler.Ih_table.set_Ih_values_to_target(target_Ih_table)
      sel = scaler.Ih_table.Ih_values != 0.0
      scaler.Ih_table = scaler.Ih_table.select(sel)
      scaler.apply_selection_to_SFs(scaler.Ih_table.nonzero_weights)

  def update_for_minimisation(self, apm, curvatures=False):
    """Update the scale factors and Ih for the next iteration of minimisation."""
    for i, scaler in enumerate(self.unscaled_scalers):
      basis_fn = basis_function(apm.apm_list[i]).return_basis()
      apm.apm_list[i].derivatives = basis_fn[1]
      scaler.Ih_table.inverse_scale_factors = basis_fn[0]

  def join_multiple_datasets(self):
    '''method to create a joint reflection table'''
    scalers = []
    if not self.params.scaling_options.target_intensities:
      scalers.extend(self.single_scalers)
    scalers.extend(self.unscaled_scalers)
    super(TargetScaler, self).join_datasets_from_scalers(scalers)

  def perform_scaling(self):
    super(TargetScaler, self).perform_scaling(target_type=ScalingTargetFixedIH)


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
    self._Ih_table = SingleIhTable(self._reflection_table, self.space_group, weights)
    logger.info('NullScaler contains %s reflections', n_refl)
    logger.info('Completed configuration of NullScaler. \n\n' + '='*80 + '\n')

  def expand_scales_to_all_reflections(self, caller=None, calc_cov=False):
    pass

  def perform_error_optimisation(self):
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
