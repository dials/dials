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
  some bookkeeping to allow selection of the parameters for one component.
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
  for multiple datasets that are simultaneously being minimised.'''
  '''Initialise with two lists of components and selections, each item of which
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
