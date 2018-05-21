"""
Tests for the scaling restraints module.
"""
from collections import OrderedDict
import pytest
from mock import Mock
from scitbx import sparse
from dials.array_family import flex
from dials.algorithms.scaling.scaling_restraints import \
  ScalingRestraints, MultiScalingRestraints

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
def mock_parameter_manager(mock_restrained_component,
  mock_unrestrained_component):
  """Mock a parameter manager to handle the components required for the
  ScalingRestraints class."""
  apm = Mock()
  apm.components = OrderedDict({'restrained' :
    {'object': mock_restrained_component,
    'n_params': mock_restrained_component.n_params, 'start_idx' : 0}})
  apm.components.update({'unrestrained' :
    {'object' : mock_unrestrained_component,
    'n_params': mock_unrestrained_component.n_params,
    'start_idx' : mock_restrained_component.n_params}})
  apm.n_active_params = (mock_restrained_component.n_params +
    mock_unrestrained_component.n_params)
  return apm

@pytest.fixture
def mock_unrestrained_apm(mock_unrestrained_component):
  """Mock a paramater manager to handle no restrained components."""
  apm = Mock()
  apm.components = OrderedDict({'unrestrained' :
    {'object' : mock_unrestrained_component,
    'n_params': mock_unrestrained_component.n_params,
    'start_idx' : 0}})
  apm.n_active_params = mock_unrestrained_component.n_params
  return apm

@pytest.fixture
def mock_multi_apm(mock_parameter_manager):
  """Mock a multi-dataset parameter manager."""
  multi_apm = Mock()
  multi_apm.apm_list = [mock_parameter_manager, mock_parameter_manager]
  n_params = mock_parameter_manager.n_active_params
  multi_apm.n_active_params = n_params * 2
  multi_apm.apm_data = {0 : {'start_idx' : 0, 'end_idx' : n_params},
    1 : {'start_idx': n_params, 'end_idx' : n_params * 2}}
  return multi_apm

@pytest.fixture
def mock_multi_unrestrained_apm(mock_unrestrained_apm):
  """Mock a paramater manager to handle no restrained components."""
  multi_apm = Mock()
  multi_apm.apm_list = [mock_unrestrained_apm, mock_unrestrained_apm]
  return multi_apm

def test_unrestrained_ScalingRestraints(mock_unrestrained_apm,
    mock_multi_unrestrained_apm):
  """Test the case of unrestrained components. None should be returned in each
  case."""

  assert ScalingRestraints(mock_unrestrained_apm).calculate_restraints() is None
  assert ScalingRestraints(mock_unrestrained_apm
    ).calculate_jacobian_restraints() is None

  assert MultiScalingRestraints(mock_multi_unrestrained_apm
    ).calculate_restraints() is None
  assert MultiScalingRestraints(mock_multi_unrestrained_apm
    ).calculate_jacobian_restraints() is None


def test_ScalingRestraints(mock_parameter_manager, mock_restrained_component,
    mock_unrestrained_component):
  """Test for the single scaling restraints manager."""

  # Test the call to calculate restraints. This should return a residual
  # vector of the same length as the restraints of the restrained component,
  # and a gradient vector of the total length of all parameters.
  restraints = ScalingRestraints(mock_parameter_manager).calculate_restraints()
  abs_restraints = mock_restrained_component.calculate_restraints()
  assert list(restraints[0]) == list(abs_restraints[0])
  assert list(restraints[1]) == (list(abs_restraints[1]) +
    [0.0] * mock_unrestrained_component.n_params)

  # Test the call to calculate jacobian restraints. This should return a
  # restraints vector, jacobian and weights vector from the components. The
  # restraints and weights should be the length of the restrained component.
  # The jacobian has n_rows equal to the number of restrainted parameters,
  # n_cols equal to the total number of parameters. Check that these are
  # correctly composed.
  jacobian_restraints = ScalingRestraints(mock_parameter_manager
    ).calculate_jacobian_restraints()
  abs_restraints = mock_restrained_component.calculate_jacobian_restraints()
  assert list(jacobian_restraints[0]) == list(abs_restraints[0])
  assert jacobian_restraints[1].n_cols == mock_parameter_manager.n_active_params
  assert jacobian_restraints[1].n_rows == mock_restrained_component.n_params
  for i in range(mock_restrained_component.n_params):
    for j in range(mock_restrained_component.n_params):
      assert jacobian_restraints[1][i, j] == abs_restraints[1][i, j]
  # All other elements should be zero
  assert abs_restraints[1].non_zeroes == jacobian_restraints[1].non_zeroes

def test_MultiScalingRestraints(mock_multi_apm, mock_restrained_component,
    mock_unrestrained_component):
  """Test for the multi-dataset scaling restraints manager."""

  # Test the call to calculate restraints. Expected return is the individual
  # dataset vectors joined together.
  restraints = MultiScalingRestraints(mock_multi_apm).calculate_restraints()
  abs_restraints = mock_restrained_component.calculate_restraints()
  assert list(restraints[0]) == (list(abs_restraints[0]) +
    list(abs_restraints[0]))
  assert list(restraints[1]) == (list(abs_restraints[1]) +
    [0.0] * mock_unrestrained_component.n_params + list(abs_restraints[1]) +
    [0.0] * mock_unrestrained_component.n_params)

  # Test the call to calculate jacobian restraints. Again, the expected return
  # is the individual dataset vectors joined together.
  jacobian_restraints = MultiScalingRestraints(mock_multi_apm
    ).calculate_jacobian_restraints()
  abs_restraints = mock_restrained_component.calculate_jacobian_restraints()

  assert list(jacobian_restraints[0]) == (list(abs_restraints[0]) +
    list(abs_restraints[0]))
  n_abs_params = mock_restrained_component.n_params
  n_total_params = mock_multi_apm.apm_data[0]['end_idx']
  assert jacobian_restraints[1].n_cols == mock_multi_apm.n_active_params
  assert jacobian_restraints[1].n_rows == n_abs_params * 2
  # Check that both restraints jacobians were set in correct location.
  for i in range(mock_restrained_component.n_params):
    for j in range(mock_restrained_component.n_params):
      assert jacobian_restraints[1][i, j] == abs_restraints[1][i, j]
      assert jacobian_restraints[1][i + n_abs_params, j + n_total_params] == (
        abs_restraints[1][i, j])
  assert abs_restraints[1].non_zeroes * 2 == jacobian_restraints[1].non_zeroes
