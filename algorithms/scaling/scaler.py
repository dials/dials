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
import time
from dials.array_family import flex
from cctbx import crystal, sgtbx
from scitbx import sparse
from dials_scaling_ext import row_multiply
#from libtbx import easy_mp
from libtbx.table_utils import simple_table
from dials.algorithms.scaling.basis_functions import basis_function
from dials.algorithms.scaling.scaling_utilities import calc_normE2
from dials.algorithms.scaling.outlier_rejection import reject_outliers
from dials.algorithms.scaling.Ih_table import IhTable
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
from dials_scaling_ext import calc_sigmasq as cpp_calc_sigmasq
logger = logging.getLogger('dials')


class ScalerBase(object):
  """Base class for all Scalers (single and multiple)."""

  __metaclass__ = abc.ABCMeta

  def __init__(self):
    self._experiments = None
    self._space_group = None
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
    assert new_Ih_table.id_ == 'IhTable', msg
    self._Ih_table = new_Ih_table

  @property
  def experiments(self):
    """The experiment object associated with the dataset."""
    return self._experiments

  @property
  def space_group(self):
    """The space group associated with the dataset."""
    return self._space_group

  @space_group.setter
  def space_group(self, new_sg):
    """Take in either a space_group symbol or object """
    if self._space_group:
      current_sg = self._space_group.info()
    else:
      current_sg = None
    if isinstance(new_sg, str):
      crystal_symmetry = crystal.symmetry(space_group_symbol=new_sg)
      self._space_group = crystal_symmetry.space_group()
    elif isinstance(new_sg, sgtbx.space_group):
      self._space_group = new_sg
      new_sg = new_sg.info()
    else:
      raise AssertionError('''Space group not recognised as a space group symbol
        or cctbx.sgtbx.space group object.''')
    if self.experiments:
      self.experiments.crystal.set_space_group(self._space_group)
    if current_sg:
      msg = ('WARNING: Manually overriding space group from {0} to {1}. {sep}'
      'If the reflection indexing in these space groups is different, {sep}'
      'bad things may happen!!! {sep}').format(current_sg, new_sg, sep='\n')
      logger.info(msg)

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

  @abc.abstractmethod
  def expand_scales_to_all_reflections(self, caller=None, calc_cov=False):
    """Expand scales from a subset to all reflections."""

  @classmethod
  def _scaling_subset(cls, reflection_table, params):
    """Select reflections with non-zero weight and update scale weights."""
    sel = ~(reflection_table.get_flags(reflection_table.flags.bad_for_scaling,
      all=False))
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
      st = time.time()
      refinery = scaling_refinery(engine=engine, target=target_type(self, apm),
        prediction_parameterisation=apm, max_iterations=max_iterations)
      refinery.run()
      ft = time.time()
      logger.info("Time taken for refinement %s", (ft - st))
      self = refinery.return_scaler()
      logger.info('\n'+'='*80+'\n')

  @abc.abstractmethod
  def perform_error_optimisation(self):
    """Optimise an error model."""


class SingleScalerBase(ScalerBase):
  """
  Parent class for single-dataset Scalers, containing a standard
  setup routine for the reflection_table - takes in params, experiment
  and reflection.
  """

  id_ = 'single'

  def __init__(self, params, experiment, reflection, scaled_id=0, for_multi=False):
    logger.info('Configuring a Scaler for a single dataset. \n')
    logger.info('The dataset id assigned to this reflection table is %s, \n'
      'the scaling model type being applied is %s. \n', scaled_id,
      experiment.scaling_model.id_)
    super(SingleScalerBase, self).__init__()
    self._experiments = experiment
    self._params = params
    self._space_group = self.experiments.crystal.get_space_group()
    if self._params.scaling_options.space_group:
      self.space_group = self._params.scaling_options.space_group
    n_model_params = sum([val.n_params for val in self.components.itervalues()])
    self._var_cov = sparse.matrix(n_model_params, n_model_params)
    # Rename the id for now so that we have unique dataset ids
    reflection['id'] = flex.int([scaled_id] * reflection.size())
    self._initial_keys = [key for key in reflection.keys()]
    # Choose intensities, map to asu, assign unique refl. index
    reflection_table = self._reflection_table_setup(self._initial_keys,
      reflection) #flag not integrated, bad d, bad partiality
    reflection_table = self._select_optimal_intensities(reflection_table,
      self.params.scaling_options.integration_method) #flag bad variance,
    # Calculate values for later filtering, but don't filter here!!!
    reflection_table = calc_normE2(reflection_table, self.experiments)
    #if 'centric_flag' in reflection_table: #If E2 values have been calculated
    #  reflection_table = set_wilson_outliers(reflection_table)
    if self.params.scaling_options.outlier_rejection != '0':
      reflection_table = self.round_of_outlier_rejection(reflection_table)
    self._reflection_table = reflection_table
    self._configure_reflection_table()
    self.select_reflections_for_scaling(for_multi=for_multi)
    #if self.params.scaling_options.n_proc > 1:
    #  self.Ih_table.split_into_blocks(self.params.scaling_options.n_proc)
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
    basis_fn = basis_function(apm, curvatures).calculate_scales_and_derivatives()
    if self.Ih_table.free_set_sel:
      dervis_trans = basis_fn[1][0].transpose() #what about when also blocked?
      work_set = ~self.Ih_table.free_set_sel
      derivs = dervis_trans.select_columns(work_set.iselection())
      self.Ih_table.derivatives = [derivs.transpose()]
    else:
      self.Ih_table.set_derivatives(basis_fn[1])
    if curvatures:
      apm.curvatures = basis_fn[2]
    self.Ih_table.set_inverse_scale_factors(basis_fn[0])
    self.Ih_table.update_weights()
    self.Ih_table.calc_Ih()

  @staticmethod
  def _reflection_table_setup(initial_keys, reflections):
    """Initial filter to select integrated reflections."""
    mask = ~reflections.get_flags(reflections.flags.integrated)
    d_mask = reflections['d'] <= 0.0
    partials_mask = reflections['partiality'] < 0.6
    reflections.set_flags(mask or partials_mask or d_mask,
      reflections.flags.excluded_for_scaling)
    logger.info('%s reflections not suitable for scaling (low partiality,\n'
      'not integrated etc).\n', reflections.get_flags(
      reflections.flags.excluded_for_scaling).count(True))
    if not 'inverse_scale_factor' in initial_keys:
      reflections['inverse_scale_factor'] = (flex.double(reflections.size(), 1.0))
    return reflections

  @classmethod
  def _select_optimal_intensities(cls, reflection_table, integration_method):
    """Choose which intensities to use for scaling."""
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
      if 'dqe' in reflection_table:
        conversion /= reflection_table['dqe']
      if 'qe' in reflection_table:
        conversion /= reflection_table['qe']
      reflection_table['intensity'] = (
        reflection_table['intensity.'+intstr+'.value'] * conversion)
      reflection_table['variance'] = (
        reflection_table['intensity.'+intstr+'.variance'] * (conversion**2))
      logger.info('%s intensity values will be used for scaling. \n',
        'Profile fitted' if intstr == 'prf' else 'Summation integrated')
    #perform a combined prf/sum in a similar fashion to aimless
    else:
      logger.info('Intensity selection choice does not have an implementation \n'
        'using prf values.')
      return cls._select_optimal_intensities(reflection_table, 'prf')
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
      scales = flex.double([])
      for inv_scale in component.inverse_scales:
        scales.extend(inv_scale)
      scaled_invsf *= scales
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

  def update_error_model(self, error_model):
    """Apply a correction to try to improve the error estimate."""
    self.Ih_table.update_error_model(error_model.refined_parameters)
    self.experiments.scaling_model.set_error_model(error_model)

  def compute_restraints_residuals_jacobian(self, apm):
    """Calculate a restraint for the jacobian."""
    return ScalingRestraints(apm).calculate_jacobian_restraints()

  def calculate_restraints(self, apm):
    """Calculate a restraint for the residual and gradient."""
    return ScalingRestraints(apm).calculate_restraints()

  def round_of_outlier_rejection(self, reflection_table):
    """calculate outliers from the reflections in the Ih_table,
    and use these to filter the reflection table and Ih_table."""
    reflection_table = reject_outliers(reflection_table, self.space_group,
      self.params.scaling_options.outlier_rejection,
      self.params.scaling_options.outlier_zmax)
    msg = ('A round of outlier rejection has been performed, in total {0} {sep}'
        'outliers have now been identified. {sep}'.format(
        reflection_table.get_flags(reflection_table.flags.outlier_in_scaling
        ).count(True), sep='\n'))
    logger.info(msg)
    return reflection_table

  def apply_error_model_to_variances(self):
    """Apply an aimless-like error model to the variances."""
    error_model = self.experiments.scaling_model.error_model
    self.reflection_table['variance'] = error_model.update_variances(
      self.reflection_table['variance'])
    logger.info('The error model determined has been applied to the variances')

  '''def apply_selection_to_SFs(self, sel):
    """Updates data from within current SF objects using the given selection.
       Required for targeted scaling."""
    for component in self.components.itervalues():
      component.update_reflection_data(self.reflection_table, sel)'''

  def select_reflections_for_scaling(self, for_multi=True):
    """Select a subset of reflections, create and Ih table and update the
    model components."""
    selection = self._scaling_subset(self.reflection_table, self.params)
    self.scaling_selection = selection
    if not for_multi:
      self._Ih_table = IhTable([(self.reflection_table, selection)],
        self.space_group, n_blocks=self.params.scaling_options.n_proc,
        weighting_scheme=self.params.weighting.weighting_scheme)
      block_selections = self.Ih_table.blocked_selection_list[0]
      '''if self.params.scaling_options.use_free_set and split_free_set:
        self._Ih_table.split_into_free_work(
          self.params.scaling_options.free_set_percentage)'''
      for component in self.components.itervalues():
        component.update_reflection_data(self.reflection_table, selection, block_selections)

  def perform_error_optimisation(self):
    self.Ih_table = IhTable([(self._reflection_table, None)], self.space_group)
    refinery = error_model_refinery(engine='SimpleLBFGS',
      target=ErrorModelTarget(BasicErrorModel(self.Ih_table.blocked_data_list[0])),
      max_iterations=100)
    refinery.run()
    error_model = refinery.return_error_manager()
    self.update_error_model(error_model)

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
    self._space_group = single_scalers[0].space_group
    self.active_scalers = None
    self._initial_keys = self.single_scalers[0].initial_keys
    self._params = params
    self._experiments = experiments[0]

  def calculate_restraints(self, apm):
    """Calculate a restraints residuals/gradient vector for multiple datasets."""
    return MultiScalingRestraints(apm).calculate_restraints()

  def compute_restraints_residuals_jacobian(self, apm):
    """Calculate a restraints residuals/jacobian for multiple datasets."""
    return MultiScalingRestraints(apm).calculate_jacobian_restraints()

  def expand_scales_to_all_reflections(self, caller=None, calc_cov=False):
    for scaler in self.active_scalers:
      scaler.expand_scales_to_all_reflections(caller=self, calc_cov=calc_cov)

  @abc.abstractmethod
  def join_multiple_datasets(self):
    """Combine all datasets into a single reflection table."""

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
      self.params.scaling_options.outlier_rejection,
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

  def update_for_minimisation(self, apm, curvatures=False):
    """Update the scale factors and Ih for the next iteration of minimisation."""
    scales_list = []
    derivs_list = []
    for apm_i in apm.apm_list:
      basis_fn = basis_function(apm_i, curvatures).calculate_scales_and_derivatives()
      scales_list.append(basis_fn[0])
      derivs_list.append(basis_fn[1])
    scales = zip(*scales_list)
    scale_list_for_Ih = []
    for scales_tuple in scales:
      scales_for_block = flex.double()
      for array in scales_tuple:
        scales_for_block.extend(array)
      scale_list_for_Ih.append(scales_for_block)
    derivs = zip(*derivs_list)
    derivs_list_for_Ih = []
    for i, derivs_tuple in enumerate(derivs):
      deriv_matrix = sparse.matrix(scale_list_for_Ih[i].size(), apm.n_active_params)
      start_row_no = 0
      for j, deriv in enumerate(derivs_tuple):
        deriv_matrix.assign_block(deriv, start_row_no, apm.apm_data[j]['start_idx'])
        start_row_no += deriv.n_rows
      derivs_list_for_Ih.append(deriv_matrix)
    self.Ih_table.set_inverse_scale_factors(scale_list_for_Ih)
    self.Ih_table.set_derivatives(derivs_list_for_Ih)
    self.Ih_table.update_weights()
    # The parallelisation below would work if sparse matrices were
    # pickleable (I think!) - with more benefit for larger number of datasets.
    '''def task_wrapper(block):
      bf = basis_function(block)
      s, d = bf.calculate_scales_and_derivatives()
      return s, d
    blocks = apm.apm_list
    task_results = easy_mp.parallel_map(func=task_wrapper, iterable=blocks,
      processes=n_datasets, method="multiprocessing",
      preserve_exception_message=True)
    scales_list, derivs_list = zip(*task_results)'''


  def select_reflections_for_scaling(self):
    """Select reflections for scaling in individual scalers."""
    #Update the data from the Ih_table
    for scaler in self.active_scalers:
      scaler.select_reflections_for_scaling(for_multi=True)

  def perform_error_optimisation(self):
    """Perform an optimisation of the sigma values."""
    for s_scaler in self.active_scalers:
      s_scaler.Ih_table = IhTable([(s_scaler.reflection_table, None)],
        s_scaler.space_group)
      refinery = error_model_refinery(engine='SimpleLBFGS',
        target=ErrorModelTarget(BasicErrorModel(
          s_scaler.Ih_table.blocked_data_list[0])), max_iterations=100)
      refinery.run()
      error_model = refinery.return_error_manager()
      s_scaler.update_error_model(error_model)
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
    logger.info('Determining symmetry equivalent reflections across datasets.\n')
    self._Ih_table = IhTable([(x.reflection_table, x.scaling_selection)
      for x in self.single_scalers], self._space_group,
      n_blocks=self.params.scaling_options.n_proc)
    '''if self.params.scaling_options.use_free_set:
      self._Ih_table.split_into_free_work(
        self.params.scaling_options.free_set_percentage)'''
    self.active_scalers = self.single_scalers
    #now add data to scale components from datasets
    for i, scaler in enumerate(self.active_scalers):
      for component in scaler.components.itervalues():
        component.update_reflection_data(scaler.reflection_table,
          scaler.scaling_selection, self.Ih_table.blocked_selection_list[i])
    logger.info('Completed configuration of MultiScaler. \n\n' + '='*80 + '\n')

  def update_error_model(self, error_model):
    """Update the error model in Ih table."""
    self.Ih_table.update_error_model(error_model)

  def update_for_minimisation(self, apm, curvatures=False):
    '''update the scale factors and Ih for the next iteration of minimisation,
    update the x values from the amp to the individual apms, as this is where
    basis functions, target functions etc get access to the parameters.'''
    super(MultiScaler, self).update_for_minimisation(apm, curvatures)
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
    logger.info('Determining symmetry equivalent reflections across datasets.\n')
    self.unscaled_scalers = unscaled_scalers
    self.active_scalers = self.unscaled_scalers
    self._target_Ih_table = IhTable([(x.reflection_table, x.scaling_selection)
      for x in self.single_scalers], self._space_group,
      n_blocks=1)#Keep in one table for matching below
    self._Ih_table = IhTable([(x.reflection_table, x.scaling_selection)
      for x in self.unscaled_scalers], self._space_group,
      n_blocks=params.scaling_options.n_proc)
    self._initial_keys = self.unscaled_scalers[0].initial_keys #needed for
    # scaling against calculated Is
    self.set_Ih_values_to_target()
    #block_selections_list = self.Ih_table.blocked_Ih_table.block_selection_list
    for i, scaler in enumerate(self.active_scalers):
      for component in scaler.components.itervalues():
        component.update_reflection_data(scaler.reflection_table,
          scaler.scaling_selection, self.Ih_table.blocked_selection_list[i])
    logger.info('Completed initialisation of TargetScaler. \n' + '*'*40 + '\n')

  def set_Ih_values_to_target(self):
    """Match equivalent reflections between individual unscaled datasets and
    set target Ixh values in the Ih_table for each dataset."""
    target_Ih_table = self._target_Ih_table.blocked_data_list[0]
    for i, block in enumerate(self._Ih_table.blocked_data_list):
      block.set_Ih_values_to_target(target_Ih_table)
      sel = block.Ih_values != 0.0
      block = block.select(sel)
      self._Ih_table.apply_selection_to_selection_list(i, sel)
      '''if self.params.scaling_options.use_free_set:
        scaler.Ih_table.split_into_free_work(
          self.params.scaling_options.free_set_percentage)'''

  def join_multiple_datasets(self):
    '''method to create a joint reflection table'''
    scalers = []
    if not self.params.scaling_options.target_model:
      scalers.extend(self.single_scalers)
    scalers.extend(self.unscaled_scalers)
    super(TargetScaler, self).join_datasets_from_scalers(scalers)

  def perform_scaling(self, target_type=ScalingTargetFixedIH, engine=None,
      max_iterations=None):
    super(TargetScaler, self).perform_scaling(target_type=target_type,
      engine=engine, max_iterations=max_iterations)


class NullScaler(ScalerBase):
  """
  Null scaler to allow targeted scaling against calculated intensities.
  """

  id_ = 'null'

  def __init__(self, params, experiment, reflection, scaled_id=0):
    super(NullScaler, self).__init__()
    self._experiments = experiment
    self._params = params
    self._space_group = self.experiments.crystal.get_space_group()
    if self._params.scaling_options.space_group:
      self._space_group = self._params.scaling_options.space_group
    self._scaled_id = scaled_id
    self._reflection_table = reflection
    self._initial_keys = [key for key in self._reflection_table.keys()]
    #self._reflection_table['intensity'] = self._reflection_table[
    #  'intensity.calculated.value']
    n_refl = self._reflection_table.size()
    self._reflection_table['inverse_scale_factor'] = flex.double(n_refl, 1.0)
    self._reflection_table['variance'] = flex.double(n_refl, 1.0)
    #weights = self._reflection_table['intensity.calculated.value']
    self._reflection_table.set_flags(flex.bool(n_refl, False),
      self._reflection_table.flags.excluded_for_scaling)
    #self._Ih_table = SingleIhTable(self._reflection_table, self.space_group,
    #  weighting_scheme='unity')
    self.scaling_selection = None
    logger.info('NullScaler contains %s reflections', n_refl)
    logger.info('Completed configuration of NullScaler. \n\n' + '='*80 + '\n')

  @property
  def components(self):
    """Shortcut to scaling model components."""
    return self.experiments.scaling_model.components

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
    n_refl = components[component].inverse_scales[0].size() #should all be same
  jacobian = sparse.matrix(n_refl, n_param)
  n_cumulative_param = 0
  for component in components:
    block = components[component].derivatives[0]
    n_param = components[component].n_params
    for component_2 in components:
      if component_2 != component:
        block = row_multiply(block, components[component].inverse_scales[0])
    jacobian.assign_block(block, 0, n_cumulative_param)
    n_cumulative_param += n_param
  logger.info('Calculating error estimates of inverse scale factors. \n')
  return cpp_calc_sigmasq(jacobian.transpose(), var_cov)
