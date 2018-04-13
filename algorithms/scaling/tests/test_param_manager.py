"""
Tests for the active parameter manager module.
"""
from scitbx import sparse
from dials.array_family import flex
from active_parameter_managers import (multi_active_parameter_manager,
  active_parameter_manager, ConcurrentAPMFactory, ConsecutiveAPMFactory)
import pytest

class DummyComponent(object):
  """Dummy model component object."""
  parameters = flex.double([1.0])
  n_params = len(parameters)
  var_cov_matrix = sparse.matrix(1, 1)
  parameter_esds = None

class DummyDataManager(object):
  """Dummy class to hold a components dictionary."""
  def __init__(self, components):
    self.components = components
    self.consecutive_refinement_order = None

  def set_consecutive_order(self, orderlist):
    """Set the consecutive refinement order list."""
    self.consecutive_refinement_order = orderlist

def test_general_apm():
  """Test for the general active_parameter_manage class."""
  components = {'scale' : DummyComponent(), 'decay' : DummyComponent(),
    'absorption' : DummyComponent()}

  apm = active_parameter_manager(components, ['scale', 'decay'])
  assert 'decay' in apm.components_list
  assert 'scale' in apm.components_list
  assert 'absorption' not in apm.components_list
  assert apm.n_active_params == (components['scale'].n_params
    + components['decay'].n_params)
  n_cumul = 0
  for component in apm.components:
    assert apm.components[component]['n_params'] == components[component].n_params
    assert apm.components[component]['start_idx'] == n_cumul
    assert apm.components[component]['end_idx'] == n_cumul + apm.components[component]['n_params']
    n_cumul += apm.components[component]['n_params']

  apm.set_param_vals(flex.double([2.0, 1.5]))
  assert apm.get_param_vals() == flex.double([2.0, 1.5])
  # Test params were updated in components
  assert list(components['scale'].parameters) == [2.0]
  assert list(components['decay'].parameters) == [1.5]
  # Test selection of parameters
  decay_params = apm.select_parameters('decay')
  assert len(decay_params) == 1
  assert decay_params[0] == 1.5

  # Test calculate model state uncertainties
  var_cov = flex.double([1.0, 0.5, 0.5, 2.0])
  var_cov.reshape(flex.grid(2, 2))
  apm.calculate_model_state_uncertainties(var_cov)
  assert components['scale'].var_cov_matrix[0, 0] == 1.0
  assert components['decay'].var_cov_matrix[0, 0] == 2.0

  # Test set param esds.
  apm.set_param_esds(flex.double([0.1, 0.2]))
  assert components['scale'].parameter_esds == flex.double([0.1])
  assert components['decay'].parameter_esds == flex.double([0.2])

def test_multi_apm():
  """Test for the general multi_active_parameter_manage class."""

  components_1 = {'scale' : DummyComponent(), 'decay' : DummyComponent(),
    'absorption' : DummyComponent()}
  components_2 = {'scale' : DummyComponent(), 'decay' : DummyComponent()}

  multi_apm = multi_active_parameter_manager([components_1, components_2],
    [['scale', 'decay'], ['scale']], active_parameter_manager)

  # Test correct setup of apm_list attribute.
  for apm in multi_apm.apm_list:
    assert isinstance(apm, active_parameter_manager)
  assert len(multi_apm.apm_list) == 2
  assert multi_apm.components_list == ['scale', 'decay', 'scale']
  assert multi_apm.n_active_params == 3
  assert multi_apm.apm_data[0] == {'start_idx': 0, 'end_idx': 2}
  assert multi_apm.apm_data[1] == {'start_idx': 2, 'end_idx': 3}

  # Test parameter selection.
  multi_apm.set_param_vals(flex.double([3.0, 2.5, 2.0]))
  assert multi_apm.get_param_vals() == flex.double([3.0, 2.5, 2.0])
  assert multi_apm.select_parameters(0) == flex.double([3.0, 2.5])
  assert multi_apm.select_parameters(1) == flex.double([2.0])

  # Test setting parameter esds.
  multi_apm.set_param_esds(flex.double([0.1, 0.2, 0.3]))
  assert components_1['scale'].parameter_esds == flex.double([0.1])
  assert components_1['decay'].parameter_esds == flex.double([0.2])
  assert components_2['scale'].parameter_esds == flex.double([0.3])

  # Test setting var_cov matrices for each component.
  var_cov = flex.double([1.0, 0.5, 0.5, 0.5, 2.0, 0.5, 0.5, 0.5, 3.0])
  var_cov.reshape(flex.grid(3, 3))
  multi_apm.calculate_model_state_uncertainties(var_cov)
  assert components_1['scale'].var_cov_matrix[0, 0] == 1.0
  assert components_1['decay'].var_cov_matrix[0, 0] == 2.0
  assert components_2['scale'].var_cov_matrix[0, 0] == 3.0

def test_concurrent_apm_factory():
  """Test the apm factory for concurrent refinement."""
  components_1 = {'scale' : DummyComponent(), 'decay' : DummyComponent(),
    'absorption' : DummyComponent()}
  data_manager = DummyDataManager(components_1)

  apm_factory = ConcurrentAPMFactory([data_manager],
    apm_type=active_parameter_manager)
  apm = apm_factory.make_next_apm()
  assert isinstance(apm, active_parameter_manager)
  assert 'scale' in apm.components_list
  assert 'decay' in apm.components_list
  assert 'absorption' in apm.components_list

  components_1 = {'scale' : DummyComponent(), 'decay' : DummyComponent(),
    'absorption' : DummyComponent()}
  components_2 = {'1' : DummyComponent(), '2' : DummyComponent()}
  data_manager_1 = DummyDataManager(components_1)
  data_manager_2 = DummyDataManager(components_2)

  apm_factory = ConcurrentAPMFactory([data_manager_1, data_manager_2],
    apm_type=active_parameter_manager)
  assert apm_factory.n_cycles == 1
  multi_apm = apm_factory.make_next_apm()
  assert isinstance(multi_apm, multi_active_parameter_manager)
  for apm in multi_apm.apm_list:
    assert isinstance(apm, active_parameter_manager)
  assert 'scale' in multi_apm.apm_list[0].components_list
  assert 'decay' in multi_apm.apm_list[0].components_list
  assert 'absorption' in multi_apm.apm_list[0].components_list
  assert '1' in multi_apm.apm_list[1].components_list
  assert '2' in multi_apm.apm_list[1].components_list

@pytest.mark.skip(reason="need to use a mock object for a scaler.")
def test_consecutive_apm_factory():
  """Test the apm factory for consecutive refinement."""
  components_1 = {'scale' : DummyComponent(), 'decay' : DummyComponent(),
    'absorption' : DummyComponent()}
  data_manager = DummyDataManager(components_1)
  data_manager.set_consecutive_order([['scale', 'decay'], ['absorption']])

  '''@mock.patch("datamanager")
  def mock_datamanager(mock_class):
    mock_class.return_value.experiments.scaling_model.consecutive_refinement_order.return_value = [
      ['scale', 'decay'], ['absorption']]
    inst = datamanager()'''

  # Test single dataset case.
  apm_factory = ConsecutiveAPMFactory([data_manager],
    apm_type=active_parameter_manager)
  assert apm_factory.n_cycles == 2
  apm = apm_factory.make_next_apm()
  assert isinstance(apm, active_parameter_manager)
  assert 'scale' in apm.components_list
  assert 'decay' in apm.components_list
  assert 'absorption' not in apm.components_list
  apm = apm_factory.make_next_apm()
  assert isinstance(apm, active_parameter_manager)
  assert 'scale' not in apm.components_list
  assert 'decay' not in apm.components_list
  assert 'absorption' in apm.components_list

  # Test multi dataset case.
  components_2 = {'1' : DummyComponent(), '2' : DummyComponent()}
  data_manager_2 = DummyDataManager(components_2)
  #data_manager_2 = mock.Mock('components'=components_2,
  #  'experiments.scaling_model.consecutive_refinement_order' = [['1'], ['2']])
  data_manager_2.set_consecutive_order([['1'], ['2']])
  apm_factory = ConsecutiveAPMFactory([data_manager, data_manager_2],
    apm_type=active_parameter_manager)
  assert apm_factory.n_cycles == 2
  multi_apm = apm_factory.make_next_apm()
  assert isinstance(multi_apm, multi_active_parameter_manager)
  apm_1 = multi_apm.apm_list[0]
  assert 'scale' in apm_1.components_list
  assert 'decay' in apm_1.components_list
  assert 'absorption' not in apm_1.components_list
  assert multi_apm.apm_list[1].components_list == ['1']
  multi_apm = apm_factory.make_next_apm()
  assert isinstance(multi_apm, multi_active_parameter_manager)
  assert multi_apm.apm_list[0].components_list == ['absorption']
  assert multi_apm.apm_list[1].components_list == ['2']

  # Test multi dataset case with different number of cycles for each data_manager.
  components_2 = {'1' : DummyComponent()}
  data_manager_2 = DummyDataManager(components_2)
  data_manager_2.set_consecutive_order([['1']])
  apm_factory = ConsecutiveAPMFactory([data_manager, data_manager_2],
    apm_type=active_parameter_manager)
  assert apm_factory.n_cycles == 2
  multi_apm = apm_factory.make_next_apm()
  assert isinstance(multi_apm, multi_active_parameter_manager)
  apm_1 = multi_apm.apm_list[0]
  assert 'scale' in apm_1.components_list
  assert 'decay' in apm_1.components_list
  assert 'absorption' not in apm_1.components_list
  assert multi_apm.apm_list[1].components_list == ['1']
  multi_apm = apm_factory.make_next_apm()
  assert isinstance(multi_apm, multi_active_parameter_manager)
  assert multi_apm.apm_list[0].components_list == ['absorption']
  # Only change relative to previous test case.
  assert multi_apm.apm_list[1].components_list == []
