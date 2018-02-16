from dials.algorithms.scaling.active_parameter_managers import \
  active_parameter_manager, multi_active_parameter_manager

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
        return ConcurrentParameterFactory([scaler], mode='single',
          apm_type=scaling_active_parameter_manager)
      else:
        return ConsecutiveParameterFactory([scaler], mode='single',
          apm_type=scaling_active_parameter_manager)
    elif isinstance(scaler, Scaler.MultiScalerBase):
      if scaler.id_ == 'target':
        data_managers = scaler.unscaled_scalers
      elif scaler.id_ == 'multi':
        data_managers = scaler.single_scalers
      else:
        assert 0, 'unrecognised scaler id_ for scaler derived from MultiScalerBase'
      if scaler.params.scaling_options.concurrent_scaling:
        return ConcurrentParameterFactory(data_managers, mode='multi',
          apm_type=scaling_active_parameter_manager)
      else:
        return ConsecutiveParameterFactory(data_managers, mode='multi',
          apm_type=scaling_active_parameter_manager)
    else:
      assert 0, '''scaler not derived from Scaler.SingleScalerBase or
        Scaler.MultiScalerBase. An additional option must be defined in
        ActiveParameterFactory'''

class ConcurrentParameterFactory(object):
  '''factory to set up apms for concurrent scaling of all active parameters'''

  def __init__(self, data_managers, mode, apm_type):
    assert mode == 'single' or mode == 'multi', '''mode must be set as single
      or multi depending on number of datasets.'''
    self.data_managers = data_managers
    self.apm = None
    self.param_lists = []
    self.create_active_list(mode, apm_type)
    self.n_cycles = 1

  def create_active_list(self, mode, apm_type):
    '''return a list indicating the names of active parameters'''

    if mode == 'single':
      param_name = []
      for param in self.data_managers[0].corrections:
        param_name.append(str(param))
      if not param_name:
        assert 0, 'no parameters have been chosen for scaling, aborting process'
      self.param_lists = param_name
      self.apm = apm_type(self.data_managers[0].components, self.param_lists)

    elif mode == 'multi':
      for data_manager in self.data_managers:
        param_name = []
        for param in data_manager.corrections:
          param_name.append(str(param))
        if not param_name:
          assert 0, 'no parameters have been chosen for scaling, aborting process'
        self.param_lists.append(param_name)
      components = [i.components for i in self.data_managers]
      self.apm = multi_active_parameter_manager(components, self.param_lists,
        apm_type)

  def make_next_apm(self):
    'method to call to return the apm'
    return self.apm

class ConsecutiveParameterFactory(object):
  '''factory to set up apms for consecutive scaling of all active parameters'''

  def __init__(self, data_managers, mode, apm_type):
    assert mode == 'single' or mode == 'multi', '''mode must be set as single
      or multi depending on number of datasets.'''
    self.data_managers = data_managers
    self.mode = mode
    self.apm_type = apm_type
    self.param_lists = []
    self.n_cycles = None
    self.create_consecutive_list()

  def create_consecutive_list(self):
    '''return a list indicating the names of active parameters'''

    if self.mode == 'single':
      for cycle in self.data_managers[0].consecutive_scaling:
        corrlist = []
        for corr in cycle:
          if corr in self.data_managers[0].corrections:
            corrlist.append(corr)
        self.param_lists.append(corrlist)
      self.n_cycles = sum([1 for i in self.param_lists if i])

    elif self.mode == 'multi':
      for data_manager in self.data_managers:
        ind_param_list = []
        for cycle in data_manager.consecutive_scaling:
          corrlist = []
          for corr in cycle:
            if corr in data_manager.corrections:
              corrlist.append(corr)
          ind_param_list.append(corrlist)
        self.param_lists.append(ind_param_list)
      #now need to calculate the max number of cycles needed across all data_managers
      is_cycle_active = []
      for p_list in self.param_lists:
        for i, cycle in enumerate(p_list):
          if cycle:
            is_cycle_active.append(i)
      self.n_cycles = len(set(is_cycle_active))
      #now make sure all lists are same length
      max_len = max([len(i) for i in self.param_lists])
      for p_list in self.param_lists:
        for _ in range(max_len - len(p_list)):
          p_list.append([])

  def make_next_apm(self):
    '''generate a valid apm for minimisation (contains some active parameters,
    but not necessarily for all datasets)'''

    if self.mode == 'single':
      apm = self.apm_type(self.data_managers[0].components, self.param_lists[0])
      self.param_lists = self.param_lists[1:]
    elif self.mode == 'multi':
      params = []
      for i in range(0, len(self.param_lists)):
        params.append(self.param_lists[i][0])
        self.param_lists[i] = self.param_lists[i][1:] #remove first element
      components = [i.components for i in self.data_managers]
      apm = multi_active_parameter_manager(components, params, self.apm_type)
    if not apm.components_list: #no active parameters, so iterate
      apm = self.make_next_apm()
    return apm
