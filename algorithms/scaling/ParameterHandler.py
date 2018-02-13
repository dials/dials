import logging
from dials.array_family import flex

logger = logging.getLogger('dials')

class ActiveParameterFactory(object):

  def __init__(self, scaler):
    self.scaler = scaler
    self.apm = None
    self.param_lists = []
    self.create_active_list()

  def create_active_list(self):
    '''return a list indicating the names of active parameters'''
    #replace these methods with some kind of entry point mechanism to allow
    #scalers derived from new base scalers in future?
    from dials.algorithms.scaling import Scaler
    if isinstance(self.scaler, Scaler.SingleScalerBase):
      param_name = []
      for param in self.scaler.corrections:
        param_name.append(str(param))
      if not param_name:
        assert 0, 'no parameters have been chosen for scaling, aborting process'
      self.param_lists = param_name
      self.apm = active_parameter_manager(self.scaler, self.param_lists)

    elif isinstance(self.scaler, Scaler.MultiScalerBase):
      if self.scaler.id_ == 'target':
        scalers = self.scaler.unscaled_scalers
      elif self.scaler.id_ == 'multi':
        scalers = self.scaler.single_scalers
      else:
        assert 0, 'unrecognised scaler id_ for scaler derived from MultiScalerBase'

      for scaler in scalers:
        param_name = []
        for param in scaler.corrections:
          param_name.append(str(param))
        if not param_name:
          assert 0, 'no parameters have been chosen for scaling, aborting process'
        self.param_lists.append(param_name)
      if self.scaler.id_ == 'target':
        self.apm = target_active_parameter_manager(self.scaler, self.param_lists)
      elif self.scaler.id_ == 'multi':
        self.apm = multi_active_parameter_manager(self.scaler, self.param_lists)

    else:
      assert 0, '''scaler not derived from Scaler.SingleScalerBase or
        Scaler.MultiScalerBase. An additional option must be defined in ParameterHandler'''

  def return_apm(self):
    'method to call to return the apm'
    return self.apm

class ConsecutiveParameterFactory(object):
  def __init__(self, scaler):
    self.scaler = scaler
    self.param_lists = None
    self.create_consecutive_list()

  def create_consecutive_list(self):
    '''return a list indicating the names of active parameters'''
    #replace these methods with some kind of entry point mechanism to allow
    #scalers derived from new base scalers in future?
    from dials.algorithms.scaling import Scaler
    if isinstance(self.scaler, Scaler.SingleScalerBase):
      param_name = []
      corrections = self.scaler.experiments.scaling_model.configdict['corrections']
      if self.scaler.experiments.scaling_model.id_ == 'aimless':
        if 'scale' in corrections and 'decay' in corrections:
          param_name.append(['scale', 'decay'])
        elif 'scale' in corrections:
          param_name.append(['scale'])
        elif 'decay' in corrections:
          param_name.append(['decay'])
        if 'absorption' in self.scaler.experiments.scaling_model.configdict['corrections']:
          param_name.append(['absorption'])
      else:
        for param in corrections:
          param_name.append([str(param)])
      self.param_lists = param_name
    else:
      assert 0, 'method not yet implemented for consecutive scaling of multiple datasets.'

  def make_next_apm(self):
    apm = active_parameter_manager(self.scaler, self.param_lists[0])
    self.param_lists = self.param_lists[1:]
    return apm

class ParameterlistFactory(object):
  '''
  Factory to create parameter lists to pass to a minimiser.
  The methods should return a nested list
  i.e [['scale','decay','absorption']]
  or [['scale','decay'],['absorption']]
  '''
  @classmethod
  def full_active_list(cls, scaler):
    '''create a list with all params to include'''
    from dials.algorithms.scaling import Scaler
    if isinstance(scaler, Scaler.SingleScalerBase) or isinstance(
      scaler, Scaler.TargetScaler): # assumes single exp in targetscaler
      param_name = []
      for param in scaler.corrections:
        param_name.append('g_'+str(param))
      if not param_name:
        assert 0, 'no parameters have been chosen for scaling, aborting process'
    elif isinstance(scaler, Scaler.MultiScaler):
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
    for p_type, scalefactor in scaler.experiments.scaling_model.components.iteritems():
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
      ).format(''.join(str(i)+' ' for i in self.active_parameterisation)))


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
      ).format(''.join(str(i)+' ' for i in self.active_parameterisation)))

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
      ).format(''.join(str(i)+' ' for i in self.active_parameterisation)))
