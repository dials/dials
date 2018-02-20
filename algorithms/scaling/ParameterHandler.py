from dials.algorithms.scaling.active_parameter_managers import \
  active_parameter_manager, ConcurrentAPMFactory, ConsecutiveAPMFactory

class scaling_active_parameter_manager(active_parameter_manager):
  '''
  Adds an extra property to the apm to avoid a repetitive calculation during
  mimimisation cycles for scaling.
  '''
  def __init__(self, data_manager, selection_list):
    self.constant_g_values = None
    components = data_manager.components
    for component, obj in components.iteritems():
      if not component in selection_list:
        assert hasattr(obj, 'inverse_scales'), '''component object must have the
          attribute 'inverse_scales'.'''
        if self.constant_g_values is None:
          self.constant_g_values = obj.inverse_scales
        else:
          self.constant_g_values *= obj.inverse_scales
    super(scaling_active_parameter_manager, self).__init__(data_manager, selection_list)


class ActiveParameterFactory(object):
  '''class to create and return appropriate apm factory.'''

  @classmethod
  def create(cls, scaler):
    '''create the factory based on phil options'''
    #replace these methods with some kind of entry point mechanism to allow
    #scalers derived from new base scalers in future?
    #if scaler.params.scaling_options.concurrent_scaling:
    from dials.algorithms.scaling import Scaler
    if isinstance(scaler, Scaler.SingleScalerBase):
      if scaler.params.scaling_options.concurrent_scaling:
        return ConcurrentAPMFactory([scaler], mode='single',
          apm_type=scaling_active_parameter_manager)
      else:
        return ConsecutiveAPMFactory([scaler], mode='single',
          apm_type=scaling_active_parameter_manager)
    elif isinstance(scaler, Scaler.MultiScalerBase):
      if scaler.id_ == 'target':
        data_managers = scaler.unscaled_scalers
      elif scaler.id_ == 'multi':
        data_managers = scaler.single_scalers
      else:
        assert 0, 'unrecognised scaler id_ for scaler derived from MultiScalerBase'
      if scaler.params.scaling_options.concurrent_scaling:
        return ConcurrentAPMFactory(data_managers, mode='multi',
          apm_type=scaling_active_parameter_manager)
      else:
        return ConsecutiveAPMFactory(data_managers, mode='multi',
          apm_type=scaling_active_parameter_manager)
    else:
      assert 0, '''scaler not derived from Scaler.SingleScalerBase or
        Scaler.MultiScalerBase. An additional option must be defined in
        ActiveParameterFactory'''
