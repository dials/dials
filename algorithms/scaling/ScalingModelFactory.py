'''
Collection of factories for creating the scaling models.
'''
from collections import OrderedDict
from dials.array_family import flex
import dials.algorithms.scaling.scale_factor as SF

class Factory(object):
  '''
  Factory for creating Scaling models
  '''
  @classmethod
  def create(cls, params, experiments, reflections):
    '''
    create the scaling model defined by the params.
    '''
    for i, (exp, refl) in enumerate(zip(experiments, reflections)):
      model = experiments.scaling_models()[i]
      if model is not None:
        exp.scaling_model = model
      elif params.scaling_model == 'aimless':
        exp.scaling_model = AimlessSMFactory.create(params, exp, refl)
      elif params.scaling_model == 'xscale':
        exp.scaling_model = XscaleScalingModel.create(params, exp, refl)
      elif params.scaling_model == 'KB':
        exp.scaling_model = KBSMFactory.create()
      else:
        assert 0, 'scaling model not recognised'
    return experiments

class AimlessSMFactory(object):
  '''
  Factory for creating an aimless-like scaling model.
  '''
  @classmethod
  def create(cls, params, experiments, reflections):
    '''
    create the scaling model defined by the params.
    '''
    scale_rot_int = params.parameterisation.rotation_interval + 0.001
    osc_range = experiments.scan.get_oscillation_range()
    one_osc_width = experiments.scan.get_oscillation()[1]
    if ((osc_range[1] - osc_range[0])/ scale_rot_int) % 1 < 0.33:
      #for scale and decay, if last bin less than 33% filled, increase rot_int
      #and extend by 0.001 to make sure all datapoints within min/max''
      n_phi_bins = int((osc_range[1] - osc_range[0])/ scale_rot_int)
      scale_rot_int = (osc_range[1] - osc_range[0])/float(n_phi_bins) + 0.001
    s_norm_fac = 0.9999 * one_osc_width / scale_rot_int
    na = reflections['xyzobs.px.value'].parts()[2] * s_norm_fac
    #need one parameter more extremal than the max/min norm values at each side
    n_scale_param = int(max(na)//1) - int(min(na)//1) + 3
    scale_parameters = flex.double([1.0] * n_scale_param)

    decay_rot_int = params.parameterisation.B_factor_interval + 0.001
    if ((osc_range[1] - osc_range[0]) / decay_rot_int) % 1 < 0.33:
      n_phi_bins = int((osc_range[1] - osc_range[0]) / decay_rot_int)
      decay_rot_int = (osc_range[1] - osc_range[0])/float(n_phi_bins) + 0.001
    d_norm_fac = 0.9999 * one_osc_width / decay_rot_int
    nt = reflections['xyzobs.px.value'].parts()[2] * d_norm_fac
    n_decay_param = int(max(nt)//1) - int(min(nt)//1) + 3
    decay_parameters = flex.double([0.0] * n_decay_param)

    lmax = params.parameterisation.lmax
    n_abs_param = (2*lmax) + (lmax**2)  #arithmetic sum formula (a1=3, d=2)
    abs_parameters = flex.double([0.0] * n_abs_param)

    return AimlessScalingModel(s_norm_fac, scale_parameters, d_norm_fac,
      decay_parameters, abs_parameters)

class XscaleSMFactory(object):
  '''
  Factory for creating an scale-like scaling model.
  '''
  @classmethod
  def create(cls, params, experiments, reflections):
    '''create an XScale scaling model.'''
    assert 0, 'method not yet implemented'

class KBSMFactory(object):
  '''
  Factory for creating a KB scaling model.
  '''
  @classmethod
  def create(cls):
    '''create the simple KB scaling model.'''
    K_param = flex.double([1.0])
    B_param = flex.double([1.0])
    return KBScalingModel(K_param, B_param)


class ScalingModelBase(object):
  '''Base class for Scale Factories'''
  def __init__(self):
    self._scaling_model = None
    self._components = OrderedDict()
    self._id_ = None
  
  @property
  def id_(self):
    return self._id_

  @property
  def components(self):
    'components of the model, a dictionary'
    return self._components


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
