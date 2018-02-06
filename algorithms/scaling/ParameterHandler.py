import logging
from dials.array_family import flex

logger = logging.getLogger('dials')

class ActiveParameterFactory(object):

  def __init__(self, scaler):
    self.scaler = scaler
    self.param_lists = []
    self.create_active_list()

  def create_active_list(self):
    from dials.algorithms.scaling import ScalerFactory
    if isinstance(self.scaler, ScalerFactory.SingleScaler):
      param_name = []
      for param in self.scaler.corrections:
        param_name.append('g_'+str(param))
      if not param_name:
        assert 0, 'no parameters have been chosen for scaling, aborting process'
      self.param_lists.append(param_name)

    if isinstance(self.scaler, ScalerFactory.TargetScaler):
      for scaler in self.scaler.unscaled_scalers:
        param_name = []
        for param in scaler.corrections:
          param_name.append('g_'+str(param))
        if not param_name:
          assert 0, 'no parameters have been chosen for scaling, aborting process'
        self.param_lists.append(param_name)    
      #param_name = []
      #for param in self.scaler.dm1.corrections:# assumes single exp in targetscaler
      #  param_name.append('g_'+str(param))
      #if not param_name:
      #  assert 0, 'no parameters have been chosen for scaling, aborting process'
      #self.param_lists.append(param_name)

    elif isinstance(self.scaler, ScalerFactory.MultiScaler):
      for scaler in self.scaler.single_scalers:
        param_name = []
        for param in scaler.corrections:
          param_name.append('g_'+str(param))
        if not param_name:
          assert 0, 'no parameters have been chosen for scaling, aborting process'
        self.param_lists.append(param_name)

  def return_active_list(self):
    return self.param_lists




class ParameterlistFactory(object):
  '''
  Factory to create parameter lists to pass to a minimiser.
  The methods should return a nested list
  i.e [['g_scale','g_decay','g_absorption']]
  or [['g_scale','g_decay'],['g_absorption']]
  '''
  @classmethod
  def full_active_list(cls, scaler):
    '''create a list with all params to include'''
    from dials.algorithms.scaling import ScalerFactory
    if isinstance(scaler, ScalerFactory.SingleScaler) or isinstance(
      scaler, ScalerFactory.TargetScaler): # assumes single exp in targetscaler
      param_name = []
      for param in scaler.corrections:
        param_name.append('g_'+str(param))
      if not param_name:
        assert 0, 'no parameters have been chosen for scaling, aborting process'
    elif isinstance(scaler, ScalerFactory.MultiScaler):
      param_name_list = []
      for scaler in scaler.single_scalers:
        param_name = []
        for param in scaler.corrections:
          param_name.append('g_'+str(param))
        if not param_name:
          assert 0, 'no parameters have been chosen for scaling, aborting process'
        param_name_list.append(param_name)

    return [param_name]

  @classmethod
  def consecutive_list(cls, scaler):
    '''create a nested list with all params to include'''
    param_name = []
    corrections = scaler.experiments.scaling_model.configdict['corrections']
    if scaler.experiments.scaling_model.id_ == 'aimless':
      if 'scale' in corrections and 'decay' in corrections:
        param_name.append(['g_scale', 'g_decay'])
      elif 'scale' in corrections:
        param_name.append(['g_scale'])
      elif 'decay' in corrections:
        param_name.append(['g_decay'])
      if 'absorption' in scaler.experiments.scaling_model.configdict['corrections']:
        param_name.append(['g_absorption'])
    else:
      for param in corrections:
        param_name.append(['g_'+str(param)])
    return param_name


class active_parameter_manager(object):
  ''' object to manage the current active parameters during minimisation.
  Separated out to provide a consistent interface between the scaler and
  minimiser. Takes in a scaler, needed to access SF objects through
  g_parameterisation, and a param_name list indicating the active parameters.'''
  def __init__(self, scaler, param_name):
    self.constant_g_values = None
    self.x = flex.double([])
    self.active_parameterisation = []
    self.n_params_list = [] #no of params in each SF
    self.n_cumul_params_list = [0]
    self.active_derivatives = None
    for p_type, scalefactor in scaler.g_parameterisation.iteritems():
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
  def __init__(self, multiscaler, param_list):
    self.apm_list = []
    for i, scaler in enumerate(multiscaler.single_scalers):
      self.apm_list.append(active_parameter_manager(scaler, param_list[i]))
    self.active_parameterisation = []
    for apm in self.apm_list:
      self.active_parameterisation.extend(apm.active_parameterisation)
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

class target_active_parameter_manager(object):
  ''' object to manage the current active parameters during minimisation
  for multiple datasets that are simultaneously being scaled.'''
  def __init__(self, targetscaler, param_list):
    self.apm_list = []
    for i, scaler in enumerate(targetscaler.unscaled_scalers):
      self.apm_list.append(active_parameter_manager(scaler, param_list[i]))
    self.active_parameterisation = []
    for apm in self.apm_list:
      self.active_parameterisation.extend(apm.active_parameterisation)
    self.x = flex.double([])
    self.n_params_in_each_apm = []
    self.n_cumul_params_list = [0]
    for apm in self.apm_list:
      self.x.extend(apm.x)
      self.n_params_in_each_apm.append(len(apm.x))
      self.n_cumul_params_list.append(len(self.x))
    self.n_active_params = len(self.x)
    logger.info(('Set up joint parameter handler for following corrections: {0}\n'
      ).format(''.join(i.lstrip('g_')+' ' for i in self.active_parameterisation)))