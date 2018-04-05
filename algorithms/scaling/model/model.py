from collections import OrderedDict
from dials.array_family import flex
from dials.algorithms.scaling.model.components.scale_components import \
  SingleScaleFactor, SingleBScaleFactor, SHScaleComponent
from dials.algorithms.scaling.model.components.smooth_scale_components import \
  SmoothScaleComponent1D, SmoothBScaleComponent1D, SmoothScaleComponent2D,\
  SmoothScaleComponent3D

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
        ('parameters', list(self._components[key].parameters))])})
      if self._components[key].parameter_esds:
        dictionary[key].update([('est_standard_devs',
          list(self._components[key].parameter_esds))])
    dictionary.update({'configuration_parameters' : self._configdict})
    return dictionary

  @classmethod
  def from_dict(cls, obj):
    '''create a scaling model object from a dictionary'''
    pass

  def set_error_model(self, error_model_params):
    self._configdict.update({'error_model_parameters' : error_model_params})

class PhysicalScalingModel(ScalingModelBase):
  '''Factory to create a scaling model for a physical parameterisation.'''

  id_ = 'physical'

  def __init__(self, parameters_dict, configdict, is_scaled=False):
    super(PhysicalScalingModel, self).__init__(configdict, is_scaled)
    if 'scale' in configdict['corrections']:
      scale_setup = parameters_dict['scale']
      self._components.update({'scale' : SmoothScaleComponent1D(
        scale_setup['parameters'], 'norm_rot_angle', scale_setup['parameter_esds'])})
    if 'decay' in configdict['corrections']:
      decay_setup = parameters_dict['decay']
      self._components.update({'decay' : SmoothBScaleComponent1D(
        decay_setup['parameters'], 'norm_time_values', decay_setup['parameter_esds'])})
    if 'absorption' in configdict['corrections']:
      absorption_setup = parameters_dict['absorption']
      self._components.update({'absorption' : SHScaleComponent(
        absorption_setup['parameters'], absorption_setup['parameter_esds'])})

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
    (s_params_sds, d_params_sds, a_params_sds) = (None, None, None)
    configdict = obj['configuration_parameters']
    is_scaled = obj['is_scaled']
    if 'scale' in configdict['corrections']:
      s_params = flex.double(obj['scale']['parameters'])
      if 'est_standard_devs' in obj['scale']:
        s_params_sds = flex.double(obj['scale']['est_standard_devs'])
    if 'decay' in configdict['corrections']:
      d_params = flex.double(obj['decay']['parameters'])
      if 'est_standard_devs' in obj['decay']:
        d_params_sds = flex.double(obj['decay']['est_standard_devs'])
    if 'absorption' in configdict['corrections']:
      abs_params = flex.double(obj['absorption']['parameters'])
      if 'est_standard_devs' in obj['absorption']:
        a_params_sds = flex.double(obj['absorption']['est_standard_devs'])

    parameters_dict = {
      'scale': {'parameters' : s_params, 'parameter_esds' : s_params_sds},
      'decay': {'parameters' : d_params, 'parameter_esds' : d_params_sds},
      'absorption': {'parameters' : abs_params, 'parameter_esds' : a_params_sds}}

    return cls(parameters_dict, configdict, is_scaled)


class ArrayScalingModel(ScalingModelBase):
  '''Factory to create a scaling model for an array-based parameterisation.'''

  id_ = 'array'

  def __init__(self, parameters_dict, configdict, is_scaled=False):
    super(ArrayScalingModel, self).__init__(configdict, is_scaled)
    if 'decay' in configdict['corrections']:
      decay_setup = parameters_dict['decay']
      self._components.update({'decay' : SmoothScaleComponent2D(
        decay_setup['parameters'], shape=(configdict['n_res_param'],
        configdict['n_time_param']), parameter_esds=decay_setup['parameter_esds'])})
    if 'absorption' in configdict['corrections']:
      abs_setup = parameters_dict['absorption']
      self._components.update({'absorption' : SmoothScaleComponent3D(
        abs_setup['parameters'], shape=(configdict['n_x_param'],
        configdict['n_y_param'], configdict['n_time_param']),
        parameter_esds=abs_setup['parameter_esds'])})
    if 'modulation' in configdict['corrections']:
      mod_setup = parameters_dict['modulation']
      self._components.update({'modulation' : SmoothScaleComponent2D(
        mod_setup['parameters'], shape=(configdict['n_x_mod_param'],
        configdict['n_y_mod_param']), parameter_esds=mod_setup['parameter_esds'])})

  @classmethod
  def from_dict(cls, obj):
    '''create a scaling model object from a dictionary'''
    if obj['__id__'] != cls.id_:
      raise RuntimeError('expected __id__ %s, got %s' % (cls.id_, obj['__id__']))
    configdict = obj['configuration_parameters']
    is_scaled = obj['is_scaled']
    (dec_params, abs_params, mod_params) = (None, None, None)
    (d_params_sds, a_params_sds, m_params_sds) = (None, None, None)
    if 'decay' in configdict['corrections']:
      dec_params = flex.double(obj['decay']['parameters'])
      if 'est_standard_devs' in obj['decay']:
        d_params_sds = flex.double(obj['decay']['est_standard_devs'])
    if 'absorption' in configdict['corrections']:
      abs_params = flex.double(obj['absorption']['parameters'])
      if 'est_standard_devs' in obj['absorption']:
        a_params_sds = flex.double(obj['absorption']['est_standard_devs'])
    if 'modulation' in configdict['corrections']:
      mod_params = flex.double(obj['modulation']['parameters'])
      if 'est_standard_devs' in obj['modulation']:
        m_params_sds = flex.double(obj['modulation']['est_standard_devs'])

    parameters_dict = {
      'decay': {'parameters' : dec_params, 'parameter_esds' : d_params_sds},
      'absorption': {'parameters' : abs_params, 'parameter_esds' : a_params_sds},
      'modulation': {'parameters' : mod_params, 'parameter_esds' : m_params_sds}}

    return cls(parameters_dict, configdict, is_scaled)

class KBScalingModel(ScalingModelBase):
  '''Factory to create a scaling model for a KB parameterisation.'''

  id_ = 'KB'

  def __init__(self, parameters_dict, configdict, is_scaled=False):
    super(KBScalingModel, self).__init__(configdict, is_scaled)
    if 'scale' in configdict['corrections']:
      self._components.update({'scale' : SingleScaleFactor(
        parameters_dict['scale']['parameters'],
        parameters_dict['scale']['parameter_esds'])})
    if 'decay' in configdict['corrections']:
      self._components.update({'decay' : SingleBScaleFactor(
        parameters_dict['decay']['parameters'],
        parameters_dict['decay']['parameter_esds'])})

  @classmethod
  def from_dict(cls, obj):
    '''create a scaling model object from a dictionary'''
    if obj['__id__'] != cls.id_:
      raise RuntimeError('expected __id__ %s, got %s' % (cls.id_, obj['__id__']))
    configdict = obj['configuration_parameters']
    is_scaled = obj['is_scaled']
    (s_params, d_params) = (None, None)
    (s_params_sds, d_params_sds) = (None, None)
    if 'scale' in configdict['corrections']:
      s_params = flex.double(obj['scale']['parameters'])
      if 'est_standard_devs' in obj['scale']:
        s_params_sds = flex.double(obj['scale']['est_standard_devs'])
    if 'decay' in configdict['corrections']:
      d_params = flex.double(obj['decay']['parameters'])
      if 'est_standard_devs' in obj['decay']:
        d_params_sds = flex.double(obj['decay']['est_standard_devs'])

    parameters_dict = {
      'scale': {'parameters' : s_params, 'parameter_esds' : s_params_sds},
      'decay': {'parameters' : d_params, 'parameter_esds' : d_params_sds}}

    return cls(parameters_dict, configdict, is_scaled)
