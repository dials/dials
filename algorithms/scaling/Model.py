from collections import OrderedDict
from dials.array_family import flex
import dials.algorithms.scaling.scale_factor as SF
from dxtbx.model import ScalingModelBaseIface

class ScalingModelBase(ScalingModelBaseIface):
  '''Base class for Scale Factories'''
  def __init__(self):
    #super(ScalingModelBase, self).__init__()
    self._scaling_model = None
    self._components = OrderedDict()
    self._id_ = None

  @property
  def id_(self):
    'id to store scaling model type'
    return self._id_

  @property
  def components(self):
    'components of the model, a dictionary'
    return self._components

  def to_dict(self):
    pass

  @classmethod
  def from_dict(cls, obj):
    pass


class AimlessScalingModel(ScalingModelBase):
  '''Factory to create a scaling model for an aimless-type parameterisation.'''
  def __init__(self, s_norm_fac, s_params, d_norm_fac, d_params, abs_params):
    super(AimlessScalingModel, self).__init__()
    self._id_ = 'aimless_scaling_model'
    self._components.update({'scale' : SF.SmoothScaleFactor1D(s_params)})
    self._components.update({'decay' : SF.SmoothBScaleFactor1D(d_params)})
    self._components.update({'absorption' : SF.SHScaleFactor(abs_params)})
    self._scale_normalisation_factor = s_norm_fac
    self._decay_normalisation_factor = d_norm_fac

  @property
  def scale_normalisation_factor(self):
    '''multiplicative factor to convert from pixel(z)
    to normalised rotation angle'''
    return self._scale_normalisation_factor

  @property
  def decay_normalisation_factor(self):
    '''multiplicative factor to convert from pixel(z)
    to normalised time'''
    return self._decay_normalisation_factor

  def to_dict(self):
    '''format data to dictionary for output'''
    dictionary = OrderedDict({'__id__' : self._id_})
    dictionary.update({'scale' : OrderedDict({
      'n_parameters' : len(self._components['scale'].parameters),
      'parameters' : list(self._components['scale'].parameters),
      'normalisation_factor' : self._scale_normalisation_factor})})
    dictionary.update({'decay' : OrderedDict({
      'n_parameters' : len(self._components['decay'].parameters),
      'parameters' : list(self._components['decay'].parameters),
      'normalisation_factor' : self._decay_normalisation_factor})})
    dictionary.update({'absorption' : OrderedDict({
      'n_parameters' : len(self._components['absorption'].parameters),
      'parameters' : list(self._components['absorption'].parameters)})})
    return dictionary

  @classmethod
  def from_dict(cls, obj):
    '''create a scaling model object from a dictionary'''
    if obj['__id__'] != "aimless_scaling_model":
      raise RuntimeError('expected __id__ aimless_scaling_model, got %s' % obj['__id__'])
    s_norm_fac = obj['scale']['normalisation_factor']
    s_params = flex.double(obj['scale']['parameters'])
    d_norm_fac = obj['decay']['normalisation_factor']
    d_params = flex.double(obj['decay']['parameters'])
    abs_params = flex.double(obj['absorption']['parameters'])
    return cls(s_norm_fac, s_params, d_norm_fac, d_params, abs_params)


class XscaleScalingModel(ScalingModelBase):
  '''Factory to create a scaling model for an xscale-type parameterisation.'''
  def __init__(self):
    super(XscaleScalingModel, self).__init__()
    assert 0, 'method not yet implemented'

  def to_dict(self):
    '''format data to dictionary for output'''
    assert 0, 'method not yet implemented'

  @classmethod
  def from_dict(cls, obj):
    '''create a scaling model object from a dictionary'''
    assert 0, 'method not yet implemented'


class KBScalingModel(ScalingModelBase):
  '''Factory to create a scaling model for an xscale-type parameterisation.'''
  def __init__(self, K_params, B_params):
    super(KBScalingModel, self).__init__()
    self._id_ = 'KB_scaling_model'
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
    if obj['__id__'] != "KB_scaling_model":
      raise RuntimeError('expected __id__ KB_scaling_model, got %s' % obj['__id__'])
    s_params = flex.double(obj['scale']['parameters'])
    d_params = flex.double(obj['decay']['parameters'])
    return cls(s_params, d_params)
