from collections import OrderedDict
from dials.array_family import flex
import dials.algorithms.scaling.model.scale_factor as SF

class ScalingModelBase(object):
  '''Base class for Scale Factories'''

  id_ = None

  def __init__(self, configdict, is_scaled=False):
    self._components = OrderedDict()
    self._configdict = configdict
    self._is_scaled = is_scaled

  @property
  def is_scaled(self):
    return self._is_scaled

  def set_scaling_model_as_scaled(self):
    '''update scaling model to indicate a scaling process has been performed on the data'''
    self._is_scaled = True

  def set_scaling_model_as_unscaled(self):
    '''update scaling model to show that no scaled data is associated with this model'''
    self._is_scaled = False

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
    dictionary = OrderedDict({'__id__' : self.id_})
    dictionary.update({'is_scaled' : self._is_scaled})
    for key in self.components:
      dictionary.update({key : OrderedDict([
        ('n_parameters', self._components[key].n_params),
        ('parameters', list(self._components[key].parameters)),
        ('est_standard_devs', list(self._components[key].parameter_esds))])})
    dictionary.update({'configuration_parameters' : self._configdict})
    return dictionary

  @classmethod
  def from_dict(cls, obj):
    '''create a scaling model object from a dictionary'''
    pass


class AimlessScalingModel(ScalingModelBase):
  '''Factory to create a scaling model for an aimless-type parameterisation.'''

  id_ = 'aimless'

  def __init__(self, s_params, d_params, abs_params, configdict, is_scaled=False):
    super(AimlessScalingModel, self).__init__(configdict, is_scaled)
    if 'scale' in configdict['corrections']:
      self._components.update({'scale' : SF.SmoothScaleFactor1D(s_params)})
    if 'decay' in configdict['corrections']:
      self._components.update({'decay' : SF.SmoothBScaleFactor1D(d_params)})
    if 'absorption' in configdict['corrections']:
      self._components.update({'absorption' : SF.SHScaleFactor(abs_params)})

  @property
  def scale_normalisation_factor(self):
    '''multiplicative factor to convert from pixel(z)
    to normalised rotation angle'''
    return self._configdict['s_norm_fac']

  @property
  def decay_normalisation_factor(self):
    '''multiplicative factor to convert from pixel(z)
    to normalised time'''
    return self._configdict['d_norm_fac']

  @classmethod
  def from_dict(cls, obj):
    '''create a scaling model object from a dictionary'''
    if obj['__id__'] != cls.id_:
      raise RuntimeError('expected __id__ %s, got %s' % (cls.id_, obj['__id__']))
    (s_params, d_params, abs_params) = (None, None, None)
    configdict = obj['configuration_parameters']
    is_scaled = obj['is_scaled']
    if 'scale' in configdict['corrections']:
      s_params = flex.double(obj['scale']['parameters'])
    if 'decay' in configdict['corrections']:
      d_params = flex.double(obj['decay']['parameters'])
    if 'absorption' in configdict['corrections']:
      abs_params = flex.double(obj['absorption']['parameters'])
    return cls(s_params, d_params, abs_params, configdict, is_scaled)


class XscaleScalingModel(ScalingModelBase):
  '''Factory to create a scaling model for an xscale-type parameterisation.'''

  id_ = 'xscale'

  def __init__(self, dec_params, abs_params, mod_params, configdict, is_scaled=False):
    super(XscaleScalingModel, self).__init__(configdict, is_scaled)
    if 'decay' in configdict['corrections']:
      self._components.update({'decay' : SF.SmoothScaleFactor2D(dec_params,
        shape=(configdict['n_res_param'], configdict['n_time_param']))})
    if 'absorption' in configdict['corrections']:
      self._components.update({'absorption' : SF.SmoothScaleFactor3D(abs_params,
        shape=(configdict['n_x_param'], configdict['n_y_param'],
          configdict['n_time_param']))})
    if 'modulation' in configdict['corrections']:
      self._components.update({'modulation' : SF.SmoothScaleFactor2D(mod_params,
        shape=(configdict['n_x_mod_param'], configdict['n_y_mod_param']))})

  @classmethod
  def from_dict(cls, obj):
    '''create a scaling model object from a dictionary'''
    if obj['__id__'] != cls.id_:
      raise RuntimeError('expected __id__ %s, got %s' % (cls.id_, obj['__id__']))
    configdict = obj['configuration_parameters']
    is_scaled = obj['is_scaled']
    (dec_params, abs_params, mod_params) = (None, None, None)
    if 'decay' in configdict['corrections']:
      dec_params = flex.double(obj['decay']['parameters'])
    if 'absorption' in configdict['corrections']:
      abs_params = flex.double(obj['absorption']['parameters'])
    if 'modulation' in configdict['corrections']:
      mod_params = flex.double(obj['modulation']['parameters'])
    return cls(dec_params, abs_params, mod_params, configdict, is_scaled)

class KBScalingModel(ScalingModelBase):
  '''Factory to create a scaling model for an xscale-type parameterisation.'''

  id_ = 'KB'

  def __init__(self, K_params, B_params, configdict, is_scaled=False):
    super(KBScalingModel, self).__init__(configdict, is_scaled)
    if 'scale' in configdict['corrections']:
      self._components.update({'scale' : SF.KScaleFactor(K_params)})
    if 'decay' in configdict['corrections']:
      self._components.update({'decay' : SF.BScaleFactor(B_params)})

  @classmethod
  def from_dict(cls, obj):
    '''create a scaling model object from a dictionary'''
    if obj['__id__'] != cls.id_:
      raise RuntimeError('expected __id__ %s, got %s' % (cls.id_, obj['__id__']))
    configdict = obj['configuration_parameters']
    is_scaled = obj['is_scaled']
    (s_params, d_params) = (None, None)
    if 'scale' in configdict['corrections']:
      s_params = flex.double(obj['scale']['parameters'])
    if 'decay' in configdict['corrections']:
      d_params = flex.double(obj['decay']['parameters'])
    return cls(s_params, d_params, configdict, is_scaled)
