import logging
from dials.array_family import flex

logger = logging.getLogger('dials.scale')

class active_parameter_manager(object):
  ''' object to manage the current active parameters during minimisation.
  Separated out to provide a consistent interface between the data manager and
  minimiser. Takes in a data manager, needed to access SF objects through
  g_parameterisation, and a param_name list indicating the active parameters.'''
  def __init__(self, Data_Manager, param_name):
    self.constant_g_values = None
    self.x = flex.double([])
    self.active_parameterisation = []
    self.n_params_list = [] #no of params in each SF
    self.n_cumul_params_list = [0]
    self.active_derivatives = None
    for p_type, scalefactor in Data_Manager.g_parameterisation.iteritems():
      if p_type in param_name:
        self.x.extend(scalefactor.parameters)
        self.active_parameterisation.append(p_type)
        self.n_params_list.append(scalefactor.n_params)
        self.n_cumul_params_list.append(len(self.x))
      else:
        if self.constant_g_values is None:
          self.constant_g_values = scalefactor.inverse_scales
        else:
          self.constant_g_values *= scalefactor.inverse_scales
    self.n_active_params = len(self.x)
    logger.info(('Set up parameter handler for following corrections: {0}\n'
      ).format(''.join(i.lstrip('g_')+' ' for i in self.active_parameterisation)))


class multi_active_parameter_manager(object):
  ''' object to manage the current active parameters during minimisation
  for multiple datasets that are simultaneously being scaled.'''
  def __init__(self, Data_Manager, param_name):
    self.apm_list = []
    for DM in Data_Manager.data_managers:
      self.apm_list.append(active_parameter_manager(DM, param_name))
    self.active_parameterisation = self.apm_list[0].active_parameterisation
    self.x = flex.double([])
    self.n_params_in_each_apm = []
    self.n_cumul_params_list = [0]
    self.active_derivatives = None
    for apm in self.apm_list:
      self.x.extend(apm.x)
      self.n_params_in_each_apm.append(len(apm.x))
      self.n_cumul_params_list.append(len(self.x))
    self.n_active_params = len(self.x)
    logger.info(('Set up joint parameter handler for following corrections: {0}\n'
      ).format(''.join(i.lstrip('g_')+' ' for i in self.active_parameterisation)))
