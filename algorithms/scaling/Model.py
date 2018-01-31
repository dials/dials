from collections import OrderedDict
from dials.array_family import flex
import dials.algorithms.scaling.scale_factor as SF

class ScalingModelBase(object):
  '''Base class for Scale Factories'''
  def __init__(self):
    self._components = OrderedDict()
    self._id_ = None
    self._configdict = None

  @property
  def id_(self):
    'id to store scaling model type'
    return self._id_

  @property
  def configdict(self):
    '''dictionary of configuration parameters'''
    return self._configdict

  @property
  def components(self):
    'components of the model, a dictionary'
    return self._components

  def to_dict(self):
    '''format data to dictionary for output'''
    pass

  @classmethod
  def from_dict(cls, obj):
    '''create a scaling model object from a dictionary'''
    pass


class AimlessScalingModel(ScalingModelBase):
  '''Factory to create a scaling model for an aimless-type parameterisation.'''
  def __init__(self, s_params, d_params, abs_params, configdict):
    super(AimlessScalingModel, self).__init__()
    self._id_ = 'aimless'
    self._configdict = configdict
    self._components.update({'scale' : SF.SmoothScaleFactor1D(s_params)})
    self._components.update({'decay' : SF.SmoothBScaleFactor1D(d_params)})
    self._components.update({'absorption' : SF.SHScaleFactor(abs_params)})
    #self._scale_normalisation_factor = configdict['s_norm_fac']
    #self._decay_normalisation_factor = configdict['d_norm_fac']

  @property
  def scale_normalisation_factor(self):
    '''multiplicative factor to convert from pixel(z)
    to normalised rotation angle'''
    return self._configdict['s_norm_fac']#_scale_normalisation_factor

  @property
  def decay_normalisation_factor(self):
    '''multiplicative factor to convert from pixel(z)
    to normalised time'''
    return self._configdict['d_norm_fac']#_decay_normalisation_factor

  def to_dict(self):
    '''format data to dictionary for output'''
    dictionary = OrderedDict({'__id__' : self._id_})
    dictionary.update({'scale' : OrderedDict({
      'n_parameters' : len(self.components['scale'].parameters),
      'parameters' : list(self.components['scale'].parameters)})})
    dictionary.update({'decay' : OrderedDict({
      'n_parameters' : len(self.components['decay'].parameters),
      'parameters' : list(self.components['decay'].parameters)})})
    dictionary.update({'absorption' : OrderedDict({
      'n_parameters' : len(self.components['absorption'].parameters),
      'parameters' : list(self.components['absorption'].parameters)})})
    dictionary.update({'configuration_parameters' : self._configdict})
    return dictionary

  @classmethod
  def from_dict(cls, obj):
    '''create a scaling model object from a dictionary'''
    if obj['__id__'] != "aimless":
      raise RuntimeError('expected __id__ aimless, got %s' % obj['__id__'])
    configdict = obj['configuration_parameters']
    #configdict['s_norm_fac'] = obj['scale']['normalisation_factor']
    s_params = flex.double(obj['scale']['parameters'])
    #configdict['d_norm_fac'] = obj['decay']['normalisation_factor']
    d_params = flex.double(obj['decay']['parameters'])
    abs_params = flex.double(obj['absorption']['parameters'])
    return cls(s_params, d_params, abs_params, configdict)


class XscaleScalingModel(ScalingModelBase):
  '''Factory to create a scaling model for an xscale-type parameterisation.'''
  def __init__(self, dec_params, abs_params, mod_params, configdict):
    super(XscaleScalingModel, self).__init__()
    #assert 0, "model not yet implemented"
    self._id_ = 'xscale'
    self._configdict = configdict
    self._components.update({'decay' : SF.SmoothScaleFactor2D(dec_params,
      shape=(configdict['n_time_param'], configdict['n_res_param']))})
    self._components.update({'absorption' : SF.SmoothScaleFactor2D(abs_params,
      shape=(configdict['n_time_param'], configdict['n_x_param']))})
    #self._components.update({'absorption' : SF.SmoothScaleFactor3D(abs_params,
    #  shape=(configdict['n_time_param'], configdict['n_x_param'],
    #  configdict['n_y_param']))})
    self._components.update({'modulation' : SF.SmoothScaleFactor2D(mod_params,
      shape=(configdict['n_x_mod_param'], configdict['n_y_mod_param']))})

  def to_dict(self):
    '''format data to dictionary for output'''
    dictionary = OrderedDict({'__id__' : self._id_})
    dictionary.update({'decay' : OrderedDict({
      'n_parameters' : len(self._components['decay'].parameters),
      'parameters' : list(self._components['decay'].parameters)})})
    dictionary.update({'absorption' : OrderedDict({
      'n_parameters' : len(self._components['absorption'].parameters),
      'parameters' : list(self._components['absorption'].parameters)})})
    dictionary.update({'modulation' : OrderedDict({
      'n_parameters' : len(self._components['modulation'].parameters),
      'parameters' : list(self._components['modulation'].parameters)})})
    dictionary.update({'configuration_parameters' : self._configdict})
    return dictionary

  @classmethod
  def from_dict(cls, obj):
    '''create a scaling model object from a dictionary'''
    if obj['__id__'] != "xscale":
      raise RuntimeError('expected __id__ xscale, got %s' % obj['__id__'])
    configdict = obj['configuration_parameters']
    mod_params = flex.double(obj['modulation']['parameters'])
    dec_params = flex.double(obj['decay']['parameters'])
    abs_params = flex.double(obj['absorption']['parameters'])
    return cls(dec_params, abs_params, mod_params, configdict)

class KBScalingModel(ScalingModelBase):
  '''Factory to create a scaling model for an xscale-type parameterisation.'''
  def __init__(self, K_params, B_params):
    super(KBScalingModel, self).__init__()
    self._id_ = 'KB'
    self._components.update({'scale' : SF.KScaleFactor(K_params)})
    self._components.update({'decay' : SF.BScaleFactor(B_params)})

  def to_dict(self):
    '''format data to dictionary for output'''
    dictionary = OrderedDict({'__id__' : self._id_})
    dictionary.update({'scale' : OrderedDict({
      'n_parameters' : len(self._components['scale'].parameters),
      'parameters' : list(self._components['scale'].parameters)})})
    dictionary.update({'decay' : OrderedDict({
      'n_parameters' : len(self._components['decay'].parameters),
      'parameters' : list(self._components['decay'].parameters)})})
    return dictionary

  @classmethod
  def from_dict(cls, obj):
    '''create a scaling model object from a dictionary'''
    if obj['__id__'] != "KB":
      raise RuntimeError('expected __id__ KB, got %s' % obj['__id__'])
    s_params = flex.double(obj['scale']['parameters'])
    d_params = flex.double(obj['decay']['parameters'])
    return cls(s_params, d_params)
