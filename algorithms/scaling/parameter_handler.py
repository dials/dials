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
          self.constant_g_values *= obj.inverse_scales
    super(scaling_active_parameter_manager, self).__init__(components,
      selection_list)
    self.n_obs = self.components[selection_list[0]]['object'].n_refl

def create_apm(scaler):
  """Create and return the appropriate apm factory for the scaler."""
  if scaler.id_ == 'single':
    if scaler.params.scaling_options.concurrent:
      return ConcurrentAPMFactory([scaler], scaling_active_parameter_manager)
    return ConsecutiveAPMFactory([scaler], scaling_active_parameter_manager)
  else:
    if scaler.id_ == 'target':
      data_managers = scaler.unscaled_scalers
    elif scaler.id_ == 'multi':
      data_managers = scaler.single_scalers
    else:
      assert 0, 'unrecognised scaler id_ for non-single scaler.'
    if scaler.params.scaling_options.concurrent:
      return ConcurrentAPMFactory(data_managers,
        scaling_active_parameter_manager, multi_mode=True)
    return ConsecutiveAPMFactory(data_managers,
      scaling_active_parameter_manager, multi_mode=True)
