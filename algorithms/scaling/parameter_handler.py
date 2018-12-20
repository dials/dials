"""
Extension to general active parameter manager for scaling and function
to use a scaler to determine the correct call to the apm factories.
"""
from scitbx.array_family import flex
from dials.algorithms.scaling.active_parameter_managers import \
  active_parameter_manager, ConcurrentAPMFactory, ConsecutiveAPMFactory

class scaling_active_parameter_manager(active_parameter_manager):
  """
  Adds an extra property to the apm to avoid a repetitive calculation during
  mimimisation cycles for scaling.
  """
  def __init__(self, components, selection_list):
    self.constant_g_values = None
    for component, obj in components.iteritems():
      if not component in selection_list:
        n_blocks = len(obj.n_refl)
        if self.constant_g_values is None:
          self.constant_g_values = [None] * n_blocks
          for n in range(n_blocks):
            self.constant_g_values[n] = obj.calculate_scales(n)
        else:
          for i in range(n_blocks):
            self.constant_g_values[n] *= obj.calculate_scales(n)
    super(scaling_active_parameter_manager, self).__init__(components,
      selection_list)
    n_obs = []
    for component in components:
      obs_in_component = []
      for n_refl in components[component].n_refl:
        obs_in_component.append(n_refl)
      n_obs.append(obs_in_component)
    assert all([i == n_obs[0] for i in n_obs])
    n_obs = []
    for component in components:
      n_obs.append(components[component].n_refl)
    self.n_obs = n_obs[0] #list of length n_blocks

def create_apm_factory(scaler):
  """Create and return the appropriate apm factory for the scaler.

  Supported cases - single dataset and multi/target datasets - concurrent or
  consecutive scaling, where the consecutive order is defined by
  scaler.consecutive_refinement_order."""
  if scaler.id_ == 'single':
    data_managers = [scaler]
  elif scaler.id_ == 'target' or scaler.id_ == 'multi':
    data_managers = scaler.active_scalers
  else:
    assert 0, 'unrecognised scaler id_'
  if scaler.params.scaling_options.concurrent:
    return ConcurrentAPMFactory(data_managers, scaling_active_parameter_manager)
  return ConsecutiveAPMFactory(data_managers, scaling_active_parameter_manager)
