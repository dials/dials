import logging
from collections import OrderedDict
from dials.array_family import flex

logger = logging.getLogger('dials')

class active_parameter_manager(object):
  '''
  Object to manage the current active parameters during minimisation.
  the parameter manager is initialised with components - a dict of ('name' : obj)
  pairs, and a selection list, which is a list of 'name' values to be minimised
  in this cycle. obj is the component of a model, and this code requires that
  obj has the attribute 'parameters'.
  This class stores a reference to the components to be minimised, alongside
  some bookkeeping to allow selection of the parameters for an individual component.
  '''
  def __init__(self, components, selection_list):
    self.x = flex.double([])
    self.components = OrderedDict()
    self.derivatives = None
    self.components_list = [] #just a list of the component names

    n_cumul_params = 0
    for component, obj in components.iteritems():
      if component in selection_list:
        assert hasattr(obj, 'parameters'), '''component object must have the
          attribute 'parameters' for access to the component parameters.'''
        self.x.extend(obj.parameters)
        n_params = len(obj.parameters)
        self.components.update({component : {'object' : obj, 'n_params' : n_params,
          'start_idx' : n_cumul_params, 'end_idx' : len(self.x)}})
        n_cumul_params += n_params
    self.n_active_params = len(self.x)

    for comp in self.components:
      self.components_list.extend([comp])

    logger.info(('Set up parameter manager for following corrections: {0}\n'
      ).format(''.join(str(i)+' ' for i in self.components_list)))

  def select_parameters(self, component):
    '''selects the subset of self.x corresponding to the component'''
    start_idx = self.components[component]['start_idx']
    end_idx = self.components[component]['end_idx']
    return self.x[start_idx : end_idx]

class multi_active_parameter_manager(object):
  '''
  Parameter manager to manage the current active parameters during minimisation
  for multiple datasets that are simultaneously being minimised.
  Initialise with two lists of components and selections, each item of which
  is used to initialise an active parameter manager of type apm_class.
  '''
  def __init__(self, components_list, selection_lists, apm_class):
    self.x = flex.double([])
    self.derivatives = None
    self.components_list = [] #just a list of the component names
    self.apm_list = []
    self.apm_data = OrderedDict()

    for component, selection_list in zip(components_list, selection_lists):
      self.apm_list.append(apm_class(component, selection_list))

    n_cumul_params = 0
    for i, apm in enumerate(self.apm_list):
      self.x.extend(apm.x)
      n_params = apm.n_active_params
      self.apm_data.update({i : {'start_idx': n_cumul_params, 'end_idx': len(self.x)}})
      n_cumul_params += n_params
    self.n_active_params = len(self.x)

    for apm in self.apm_list:
      for comp in apm.components:
        self.components_list.extend([comp])

    logger.info(('Set up multi-dataset parameter manager for following corrections: {0}\n'
      ).format(''.join(str(i)+' ' for i in self.components_list)))

  def select_parameters(self, apm_number):
    'selects the subset of self.x corresponding to the apm number'
    apm_data = self.apm_data[apm_number]
    return self.x[apm_data['start_idx'] : apm_data['end_idx']]

class ConcurrentAPMFactory(object):
  '''
  Factory to correctly set up a single/multi active parameter manager for
  concurrent scaling of all model components of a general data_manager
  (or multiple data managers). This extracts the name of the components from
  each data_manager and creates a list to pass on to the apm_type specified
  (e.g a subclass of active_parameter_manager). mode=single/multi returns a
  single/mutli active parameter manager.
  '''

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
      for param in self.data_managers[0].components:
        param_name.append(str(param))
      if not param_name:
        assert 0, 'no parameters have been chosen for scaling, aborting process'
      self.param_lists = param_name
      self.apm = apm_type(self.data_managers[0].components, self.param_lists)

    elif mode == 'multi':
      for data_manager in self.data_managers:
        param_name = []
        for param in data_manager.components:
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

class ConsecutiveAPMFactory(object):
  '''
  Factory to correctly set up a nested list structure to pass to
  single/multi active parameter managers for consecutive scaling of the
  model components of a general data_manager (or multiple data managers).
  Upon calling make_next_apm, the first list element for each dataset is used to
  initialise the apm_type specified (e.g a subclass of active_parameter_manager)
  and then removed from the list structure.
  make_next_apm can be called n_cycles times.
  mode=single/multi returns a single/mutli active parameter manager.
  '''
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
      for cycle in self.data_managers[0].consecutive_scaling_order:
        corrlist = []
        for corr in cycle:
          if corr in self.data_managers[0].components:
            corrlist.append(corr)
        self.param_lists.append(corrlist)
      self.n_cycles = sum([1 for i in self.param_lists if i])

    elif self.mode == 'multi':
      for data_manager in self.data_managers:
        ind_param_list = []
        for cycle in data_manager.consecutive_scaling_order:
          corrlist = []
          for corr in cycle:
            if corr in data_manager.components:
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
