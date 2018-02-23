from dials.algorithms.scaling.active_parameter_managers import \
  active_parameter_manager, ConcurrentAPMFactory, ConsecutiveAPMFactory
from dials.algorithms.scaling.scaler import SingleScalerBase,\
  MultiScalerBase

class scaling_active_parameter_manager(active_parameter_manager):
  '''
  Adds an extra property to the apm to avoid a repetitive calculation during
  mimimisation cycles for scaling.
  '''
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
    super(scaling_active_parameter_manager, self).__init__(components, selection_list)

def create_apm(scaler):
  '''method to create and return the appropriate apm factory'''
  if isinstance(scaler, SingleScalerBase):
    if scaler.params.scaling_options.concurrent:
      return ConcurrentAPMFactory([scaler], mode='single',
        apm_type=scaling_active_parameter_manager)
    else:
      return ConsecutiveAPMFactory([scaler], mode='single',
        apm_type=scaling_active_parameter_manager)
  elif isinstance(scaler, MultiScalerBase):
    if scaler.id_ == 'target':
      data_managers = scaler.unscaled_scalers
    elif scaler.id_ == 'multi':
      data_managers = scaler.single_scalers
    else:
      assert 0, 'unrecognised scaler id_ for scaler derived from MultiScalerBase'
    if scaler.params.scaling_options.concurrent:
      return ConcurrentAPMFactory(data_managers, mode='multi',
        apm_type=scaling_active_parameter_manager)
    else:
      return ConsecutiveAPMFactory(data_managers, mode='multi',
        apm_type=scaling_active_parameter_manager)
  else:
    assert 0, '''scaler not derived from Scaler.SingleScalerBase or
      Scaler.MultiScalerBase. An additional option must be defined in
      ActiveParameterFactory'''

