"""
Extension to general active parameter manager for scaling and function
to use a scaler to determine the correct call to the apm factories.
"""

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
        assert hasattr(obj, 'inverse_scales'), '''component object must have the
          attribute 'inverse_scales'.'''
        if self.constant_g_values is None:
          self.constant_g_values = obj.inverse_scales
        else:
          for i in range(len(obj.inverse_scales)):
            self.constant_g_values[i] *= obj.inverse_scales[i]
    super(scaling_active_parameter_manager, self).__init__(components,
      selection_list)
    n_obs = []
    for component in components:
      n_blocks = 0
      for n_refl in components[component].n_refl:
        n_obs.append(n_refl)
        n_blocks += 1
    assert len(set(n_obs)) == n_blocks or len(set(n_obs)) == 1
    # Assert same no of refl set in all components, or this is true and
    # same for each dataset.
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
    if scaler.params.scaling_options.concurrent:
      return ConcurrentAPMFactory([scaler], scaling_active_parameter_manager)
    return ConsecutiveAPMFactory([scaler], scaling_active_parameter_manager)
  else:
    if scaler.id_ == 'target' or scaler.id_ == 'multi':
      data_managers = scaler.active_scalers
    else:
      assert 0, 'unrecognised scaler id_ for non-single scaler.'
    if scaler.params.scaling_options.concurrent:
      return ConcurrentAPMFactory(data_managers,
        scaling_active_parameter_manager, multi_mode=True)
    return ConsecutiveAPMFactory(data_managers,
      scaling_active_parameter_manager, multi_mode=True)
