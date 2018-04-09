from dials.array_family import flex
from parameter_handler import scaling_active_parameter_manager
from active_parameter_managers import (multi_active_parameter_manager,
  active_parameter_manager)

class dummy_obj(object):
  """Dummy model component object."""
  parameters = flex.double([1.0])
  n_params = len(parameters)

def test_general_apm():
  """Test for the general active_parameter_manage class."""
  components = {'scale' : dummy_obj(), 'decay' : dummy_obj(),
    'absorption' : dummy_obj()}

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

  # Test selection of parameters

  # Test existence of methods

  # Test calculate model state uncertainties

  # Test set param esds.

