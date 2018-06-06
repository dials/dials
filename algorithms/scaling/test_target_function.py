"""
Test for the target function module.
"""
import copy
from collections import OrderedDict
import pytest
from mock import Mock, MagicMock
from scitbx import sparse
from libtbx import phil
from libtbx.test_utils import approx_equal
from dxtbx.model.experiment_list import ExperimentList
from dxtbx.model import Crystal, Scan, Beam, Goniometer, Detector, Experiment
from dials.array_family import flex
from dials.util.options import OptionParser
from dials.algorithms.scaling.scaling_library import create_scaling_model
from dials.algorithms.scaling.scaler_factory import create_scaler
from dials.algorithms.scaling.target_function import ScalingTarget, \
  ScalingTargetFixedIH
from dials.algorithms.scaling.parameter_handler import \
  scaling_active_parameter_manager
from dials.algorithms.scaling.active_parameter_managers import \
  multi_active_parameter_manager


@pytest.fixture
def large_reflection_table():
  """Create a larger reflection table"""
  return generated_10_refl()

@pytest.fixture
def small_reflection_table():
  """Create a small reflection table"""
  return generated_refl()

def generated_10_refl():
  """Generate reflection table to test the basis and target function."""
  #these miller_idx/d_values don't make physical sense, but I didn't want to
  #have to write the tests for lots of reflections.
  reflections = flex.reflection_table()
  reflections['intensity.prf.value'] = flex.double([75.0, 10.0, 100.0, 25.0,
    50.0, 100.0, 25.0, 20.0, 300.0, 10.0])
  reflections['intensity.prf.variance'] = flex.double([50.0, 10.0, 100.0, 50.0,
    10.0, 100.0, 50.0, 10.0, 100.0, 10.0])
  reflections['miller_index'] = flex.miller_index([(1, 0, 0), (0, 0, 1),
    (1, 0, 0), (1, 0, 0), (0, 0, 1), (1, 0, 0), (0, 4, 0), (0, 0, 1),
    (1, 0, 0), (0, 4, 0)]) #don't change
  reflections['d'] = flex.double([2.0, 0.8, 2.0, 2.0, 0.8, 2.0, 2.0, 0.8, 2.0,
    1.0]) #don't change
  reflections['partiality'] = flex.double(10, 1.0)
  reflections['xyzobs.px.value'] = flex.vec3_double([(0.0, 0.0, 0.0),
    (0.0, 0.0, 5.0), (0.0, 0.0, 10.0), (0.0, 0.0, 15.0), (0.0, 0.0, 20.0),
    (0.0, 0.0, 25.0), (0.0, 0.0, 30.0), (0.0, 0.0, 35.0), (0.0, 0.0, 40.0),
    (0.0, 0.0, 59.0)])
  reflections['s1'] = flex.vec3_double([(0.0, 0.1, 1.0), (0.0, 0.1, 1.0),
    (0.0, 0.1, 20.0), (0.0, 0.1, 20.0), (0.0, 0.1, 20.0), (0.0, 0.1, 20.0),
    (0.0, 0.1, 20.0), (0.0, 0.1, 20.0), (0.0, 0.1, 20.0), (0.0, 0.1, 20.0)])
  reflections.set_flags(flex.bool(10, True), reflections.flags.integrated)
  return [reflections]

def generated_refl():
  """Generate reflection table to test the basis and target function."""
  #these miller_idx/d_values don't make physical sense, but I didn't want to
  #have to write the tests for lots of reflections.
  reflections = flex.reflection_table()
  reflections['intensity.prf.value'] = flex.double([75.0, 10.0, 100.0])
  reflections['intensity.prf.variance'] = flex.double([50.0, 10.0, 100.0])
  reflections['miller_index'] = flex.miller_index([(1, 0, 0), (0, 0, 1),
    (1, 0, 0)]) #don't change
  reflections['d'] = flex.double([2.0, 0.8, 2.0]) #don't change
  reflections['lp'] = flex.double([1.0, 1.0, 1.0])
  reflections['dqe'] = flex.double([1.0, 1.0, 1.0])
  reflections['partiality'] = flex.double([1.0, 1.0, 1.0])
  reflections['xyzobs.px.value'] = flex.vec3_double([(0.0, 0.0, 0.0),
    (0.0, 0.0, 5.0), (0.0, 0.0, 10.0)])
  reflections['s1'] = flex.vec3_double([(0.0, 0.1, 1.0), (0.0, 0.1, 1.0),
    (0.0, 0.1, 20.0)])
  reflections.set_flags(flex.bool([True, True, True]),
    reflections.flags.integrated)
  return [reflections]

def generated_single_exp():
  """Generate an experiment object."""
  experiments = ExperimentList()
  exp_dict = {"__id__" : "crystal", "real_space_a": [1.0, 0.0, 0.0],
              "real_space_b": [0.0, 1.0, 0.0], "real_space_c": [0.0, 0.0, 2.0],
              "space_group_hall_symbol": " C 2y"}
  crystal = Crystal.from_dict(exp_dict)
  scan = Scan(image_range=[0, 60], oscillation=[0.0, 1.0])
  beam = Beam(s0=(0.0, 0.0, 1.01))
  goniometer = Goniometer((1.0, 0.0, 0.0))
  detector = Detector()
  experiments.append(Experiment(beam=beam, scan=scan, goniometer=goniometer,
    detector=detector, crystal=crystal))
  return experiments

def generated_param(model='KB'):
  """Generate the scaling phil param scope."""
  phil_scope = phil.parse('''
      include scope dials.algorithms.scaling.scaling_options.phil_scope
  ''', process_includes=True)

  optionparser = OptionParser(phil=phil_scope, check_format=False)
  parameters, _ = optionparser.parse_args(args=[], quick_parse=True,
    show_diff_phil=False)
  parameters.__inject__('model', model)
  parameters.parameterisation.absorption_term = False
  return parameters

@pytest.fixture
def single_exp():
  """Create an experimentlist with a single experiment."""
  return generated_single_exp()

@pytest.fixture
def physical_param():
  """Create a physical model params object."""
  return generated_param(model='physical')

@pytest.fixture
def mock_single_Ih_table():
  """Mock Ih table to use for testing the target function."""
  Ih_table = Mock()
  Ih_table.inverse_scale_factors = flex.double([1.0, 1.0/1.1, 1.0])
  Ih_table.intensities = flex.double([10.0, 10.0, 12.0])
  Ih_table.Ih_values = flex.double([11.0, 11.0, 11.0])
  # These values should give residuals of [-1.0, 0.0, 1.0]
  Ih_table.weights = flex.double([1.0, 1.0, 1.0])
  Ih_table.size = 3
  Ih_table.derivatives = sparse.matrix(3, 1)
  Ih_table.derivatives[0, 0] = 1.0
  Ih_table.derivatives[1, 0] = 2.0
  Ih_table.derivatives[2, 0] = 3.0
  Ih_table.h_index_matrix = sparse.matrix(3, 2)
  Ih_table.h_index_matrix[0, 0] = 1
  Ih_table.h_index_matrix[1, 0] = 1
  Ih_table.h_index_matrix[2, 1] = 1
  Ih_table.h_expand_matrix = Ih_table.h_index_matrix.transpose()
  return Ih_table

@pytest.fixture()
def mock_Ih_table(mock_single_Ih_table):
  """A mock Ih table for testing a target function."""
  Ih_table = MagicMock()
  Ih_table.blocked_data_list = [mock_single_Ih_table]
  Ih_table.free_Ih_table = None
  return Ih_table

@pytest.fixture()
def mock_Ih_table_workfree(mock_single_Ih_table):
  """A mock Ih table with a free set for testing a target function."""
  Ih_table = MagicMock()
  Ih_table.blocked_data_list = [mock_single_Ih_table, mock_single_Ih_table]
  Ih_table.free_Ih_table = True
  return Ih_table

@pytest.fixture
def mock_restrained_component():
  """Mock a component with restraints."""
  component = Mock()
  component.n_params = 3
  component.calculate_restraints.return_value = [flex.double([1.0, 2.0, 3.0]),
    flex.double([0.1, 0.2, 0.3])]
  jacobian_restr = sparse.matrix(component.n_params, component.n_params)
  jacobian_restr[0, 0] = 1.0
  component.calculate_jacobian_restraints.return_value = [
    flex.double([1.0, 1.1, 1.2]), jacobian_restr, flex.double([1.0, 2.0, 3.0])]
  return component

@pytest.fixture
def mock_unrestrained_component():
  """Mock a component without restraints."""
  component = Mock()
  component.n_params = 5
  component.calculate_restraints.return_value = None
  component.calculate_jacobian_restraints.return_value = None
  return component

@pytest.fixture
def mock_parameter_manager_withrestraints(mock_restrained_component):
  """Mock a parameter manager to handle the components required for the
  ScalingRestraints class."""
  apm = Mock()
  apm.components = OrderedDict({'restrained' :
    {'object': mock_restrained_component,
    'n_params': mock_restrained_component.n_params, 'start_idx' : 0}})
  apm.n_active_params = mock_restrained_component.n_params
  return apm

@pytest.fixture
def mock_parameter_manager_withoutrestraints(mock_unrestrained_component):
  """Mock a parameter manager to handle the components required for the
  ScalingRestraints class."""
  apm = Mock()
  apm.components = OrderedDict({'unrestrained' :
    {'object' : mock_unrestrained_component,
    'n_params': mock_unrestrained_component.n_params,
    'start_idx' :0}})
  apm.n_active_params = mock_unrestrained_component.n_params
  return apm

@pytest.fixture
def mock_multi_apm_withrestraints(mock_parameter_manager_withrestraints):
  """Mock a multi-dataset parameter manager."""
  apm = mock_parameter_manager_withrestraints
  multi_apm = Mock()
  multi_apm.apm_list = [apm]
  multi_apm.n_active_params = apm.n_active_params
  multi_apm.apm_data = {0 : {'start_idx' : 0, 'end_idx' : apm.n_active_params}}
  return multi_apm

@pytest.fixture
def mock_multi_apm_withoutrestraints(mock_parameter_manager_withoutrestraints):
  """Mock a multi-dataset parameter manager."""
  apm = mock_parameter_manager_withoutrestraints
  multi_apm = Mock()
  multi_apm.apm_list = [apm]
  multi_apm.n_active_params = apm.n_active_params
  multi_apm.apm_data = {0 : {'start_idx' : 0, 'end_idx' : apm.n_active_params}}
  return multi_apm

def test_target_function(mock_single_Ih_table, mock_multi_apm_withrestraints,
  mock_multi_apm_withoutrestraints):
  """Test for the ScalingTarget class."""

  # Create a scaling target and check gradients
  target = ScalingTarget()
  apm_restr = mock_multi_apm_withrestraints
  apm_norestr = mock_multi_apm_withoutrestraints
  # Below methods needed for refinement engine calls
  r, w = target.compute_residuals(mock_single_Ih_table)
  assert r.size() == w.size()

  f, g = target.compute_functional_gradients(mock_single_Ih_table)
  assert isinstance(f, float)
  assert g.size() == 1 #Number of parameters as determined by deriv matrix cols

  r, j, w = target.compute_residuals_and_gradients(mock_single_Ih_table)
  assert r.size() == w.size()
  assert j.n_cols == 1 #Number of parameters as determined by jacob matrix.
  assert j.n_rows == r.size()

  with pytest.raises(AssertionError):
    _ = target.compute_functional_gradients_and_curvatures(mock_single_Ih_table)

  restraints = target.compute_restraints_residuals_and_gradients(apm_restr)
  assert len(restraints) == 3
  assert target.param_restraints is True

  restraints = target.compute_restraints_functional_gradients_and_curvatures(apm_restr)
  assert len(restraints) == 3

  achieved = target.achieved()
  assert isinstance(achieved, bool)

  restraints = target.compute_restraints_residuals_and_gradients(apm_norestr)
  assert restraints is None
  assert target.param_restraints is False

  target = ScalingTarget() # Need to make new instance or won't calc restr as
  # param_restraints is set to False
  assert target.param_restraints is True
  restraints = target.compute_restraints_functional_gradients_and_curvatures(
    apm_norestr)
  assert restraints is None
  assert target.param_restraints is False

def test_target_rmsd_calculation(mock_Ih_table, mock_Ih_table_workfree,
  mock_multi_apm_withrestraints, mock_multi_apm_withoutrestraints):
  """Test the RMSD calculation for various scenarios - including an Ih
  table split into a work and free set and components with/without restraints.
  """
  target = ScalingTarget()
  assert len(target.rmsd_names) == 1
  assert len(target.rmsd_units) == 1
  assert target.param_restraints is True
  # with input, expect residuals of [-1, 0, 1], weights of [1, 1, 1],
  # restraints of [1, 2, 3], so expect residual of sqrt((2+6)/3)
  rmsds = target.rmsds(mock_Ih_table, mock_multi_apm_withrestraints)
  assert len(rmsds) == 1
  assert rmsds[0] == pytest.approx((8.0/3.0)**0.5, abs=1e-6)
  assert target.param_restraints is True

  rmsds = target.rmsds(mock_Ih_table_workfree, mock_multi_apm_withrestraints)
  assert len(rmsds) == 2
  assert rmsds[0] == pytest.approx((8.0/3.0)**0.5, abs=1e-6)
  assert rmsds[1] == pytest.approx((8.0/3.0)**0.5, abs=1e-6)
  assert target.param_restraints is True
  assert len(target.rmsd_names) == 2
  assert len(target.rmsd_units) == 2

  rmsds = target.rmsds(mock_Ih_table, mock_multi_apm_withoutrestraints)
  assert len(rmsds) == 1
  assert rmsds[0] == pytest.approx((2.0/3.0)**0.5, abs=1e-6)
  assert target.param_restraints is False
  assert len(target.rmsd_names) == 1
  assert len(target.rmsd_units) == 1

def test_target_fixedIh(mock_multi_apm_withoutrestraints, mock_Ih_table):
  """Test the target function for targeted scaling (where Ih is fixed)."""

  target = ScalingTargetFixedIH()
  Ih_table = mock_Ih_table.blocked_data_list[0]
  R, _ = target.compute_residuals(Ih_table)
  expected_residuals = flex.double([-1.0, 0.0, 1.0])
  assert list(R) == list(expected_residuals)
  _, G = target.compute_functional_gradients(Ih_table)
  assert list(G) == [-44.0]
  # Add in finite difference check

  J = target.calculate_jacobian(Ih_table)
  assert J.n_cols == 1
  assert J.n_rows == 3
  assert J.non_zeroes == 3
  assert J[0, 0] == -11.0
  assert J[1, 0] == -22.0
  assert J[2, 0] == -33.0

  expected_rmsd = (expected_residuals**2 / len(expected_residuals))**0.5
  assert target._rmsds is None
  target._rmsds = []
  target.rmsds(mock_Ih_table, mock_multi_apm_withoutrestraints)
  assert target._rmsds == pytest.approx([expected_rmsd])

# For testing the targetfunction calculations using finite difference methods,
# need to initialise real instances of the scaling datastructures to allow
# variation of the parameters and updating of the linked datastructures.

def test_target_gradient_calculation_finite_difference(small_reflection_table,
  single_exp, physical_param):
  """Test the calculated gradients against a finite difference calculation."""
  (test_reflections, test_experiments, params) = (
    small_reflection_table, single_exp, physical_param)
  assert len(test_experiments) == 1
  assert len(test_reflections) == 1
  experiments = create_scaling_model(params, test_experiments, test_reflections)
  scaler = create_scaler(params, experiments, test_reflections)
  assert scaler.experiments.scaling_model.id_ == 'physical'

  # Initialise the parameters and create an apm
  scaler.components['scale'].inverse_scales = flex.double([2.0, 1.0, 2.0])
  scaler.components['decay'].inverse_scales = flex.double([1.0, 1.0, 0.4])
  apm = multi_active_parameter_manager([scaler.components],
    [['scale', 'decay']], scaling_active_parameter_manager)

  # Now do finite difference check.
  target = ScalingTarget()
  scaler.update_for_minimisation(apm)
  grad = target.calculate_gradients(scaler.Ih_table.blocked_data_list[0])
  res = target.calculate_residuals(scaler.Ih_table.blocked_data_list[0])

  assert res > 1e-8, """residual should not be zero, or the gradient test
    below will not really be working!"""

  # Now compare to finite difference
  f_d_grad = calculate_gradient_fd(target, scaler, apm)
  print(list(f_d_grad))
  print(list(grad))
  assert approx_equal(list(grad), list(f_d_grad))

  sel = f_d_grad > 1e-8
  assert sel, """assert sel has some elements, as finite difference grad should
    not all be zero, or the test will not really be working!
    (expect one to be zero for KB scaling example?)"""

def test_target_jacobian_calculation_finite_difference(physical_param,
  single_exp, large_reflection_table):
  """Test the calculated jacobian against a finite difference calculation."""
  test_params, exp, test_refl = physical_param, single_exp, large_reflection_table
  test_params.parameterisation.decay_term = False
  experiments = create_scaling_model(test_params, exp, test_refl)
  assert experiments[0].scaling_model.id_ == 'physical'
  scaler = create_scaler(test_params, experiments, test_refl)

  apm = multi_active_parameter_manager([scaler.components], [['scale']],
    scaling_active_parameter_manager)

  target = ScalingTarget()
  scaler.update_for_minimisation(apm)

  fd_jacobian = calculate_jacobian_fd(target,
    scaler, apm)
  _, jacobian, _ = target.compute_residuals_and_gradients(
    scaler.Ih_table.blocked_data_list[0])

  n_rows = jacobian.n_rows
  n_cols = jacobian.n_cols

  print(jacobian)
  print(fd_jacobian)

  for i in range(0, n_rows):
    for j in range(0, n_cols):
      assert jacobian[i, j] == pytest.approx(fd_jacobian[i, j], abs=1e-4)

def calculate_gradient_fd(target, scaler, apm):
  """Calculate gradient array with finite difference approach."""
  delta = 1.0e-6
  gradients = flex.double([0.0] * apm.n_active_params)
  Ih_table = scaler.Ih_table.blocked_data_list[0]
  #iterate over parameters, varying one at a time and calculating the gradient
  for i in range(apm.n_active_params):
    new_x = copy.copy(apm.x)
    new_x[i] -= 0.5 * delta
    apm.set_param_vals(new_x)
    scaler.update_for_minimisation(apm)
    R_low = (target.calculate_residuals(Ih_table)**2) * Ih_table.weights
    new_x[i] += delta
    apm.set_param_vals(new_x)
    scaler.update_for_minimisation(apm)
    R_upper = (target.calculate_residuals(Ih_table)**2) * Ih_table.weights
    new_x[i] -= 0.5 * delta
    apm.set_param_vals(new_x)
    scaler.update_for_minimisation(apm)
    gradients[i] = (flex.sum(R_upper) - flex.sum(R_low)) / delta
  return gradients

def calculate_jacobian_fd(target, scaler, apm, block_id=0):
  """Calculate jacobian matrix with finite difference approach."""
  delta = 1.0e-8
  jacobian = sparse.matrix(scaler.Ih_table.blocked_data_list[block_id].size,
    apm.n_active_params)
  Ih_table = scaler.Ih_table.blocked_data_list[block_id]
  #iterate over parameters, varying one at a time and calculating the residuals
  for i in range(apm.n_active_params):
    new_x = copy.copy(apm.x)
    new_x[i] -= 0.5 * delta
    apm.set_param_vals(new_x)
    scaler.update_for_minimisation(apm)
    R_low = target.calculate_residuals(Ih_table)#unweighted unsquared residual
    new_x[i] += delta
    apm.set_param_vals(new_x)
    scaler.update_for_minimisation(apm)
    R_upper = target.calculate_residuals(Ih_table) #unweighted unsquared residual
    new_x[i] -= 0.5 * delta
    apm.set_param_vals(new_x)
    scaler.update_for_minimisation(apm)
    fin_difference = (R_upper - R_low) / delta
    for j in range(fin_difference.size()):
      jacobian[j, i] = fin_difference[j]
  return jacobian
