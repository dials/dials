"""
Definitions of scaling models - collections of scale components with appropriate
methods to define how these are composed into one model.
"""
import abc
import logging
import numpy as np
from collections import OrderedDict
from dials.array_family import flex
from dials.algorithms.scaling.model.components.scale_components import \
  SingleScaleFactor, SingleBScaleFactor, SHScaleComponent
from dials.algorithms.scaling.model.components.smooth_scale_components import \
  SmoothScaleComponent1D, SmoothBScaleComponent1D, SmoothScaleComponent2D,\
  SmoothScaleComponent3D

logger = logging.getLogger('dials')

class ScalingModelBase(object):
  '''Base class for Scale Factories'''

  id_ = None

  __metaclass__ = abc.ABCMeta

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

  def configure_reflection_table(self, reflection_table, experiment, params):
    return reflection_table

  def normalise_components(self):
    """Optionally define a normalisation of the parameters after scaling."""
    pass

  @property
  def configdict(self):
    '''dictionary of configuration parameters'''
    return self._configdict

  @property
  def components(self):
    'components of the model, a dictionary'
    return self._components

  @abc.abstractproperty
  def consecutive_scaling_order(self):
    """Return a nested list of correction names, to indicate the order
    to perform scaling in consecutive scaling mode if concurrent=0.
    e.g. [['scale', 'decay'], ['absorption']] would cause the first cycle to
    refine scale and decay, and then absorption in a subsequent cycle."""
    pass

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

  @property
  def consecutive_scaling_order(self):
    return [['scale', 'decay'], ['absorption']]

  def configure_reflection_table(self, reflection_table, experiment, params):
    reflection_table['phi'] = (reflection_table['xyzobs.px.value'].parts()[2]
      * experiment.scan.get_oscillation()[1])
    for component in self.components.itervalues():
      reflection_table = component.configure_reflection_table(
        reflection_table, experiment, params)
    return reflection_table

  def normalise_components(self):
    if 'scale' in self.components:
      # Do an invariant rescale of the scale at t=0 to one.'''
      sel = (self.components['scale'].normalised_values ==
        min(self.components['scale'].normalised_values))
      initial_scale = self.components['scale'].inverse_scales.select(sel)[0]
      self.components['scale'].parameters /= initial_scale
      self.components['scale'].calculate_scales_and_derivatives()
      logger.info('\nThe "scale" model component has been rescaled, so that the\n'
        'initial scale is 1.0.')
    if 'decay' in self.components:
      # Do an invariant rescale of the max B to zero.'''
      maxB = max(flex.double(np.log(self.components['decay'].inverse_scales))
                  * 2.0 * (self.components['decay'].d_values**2))
      self.components['decay'].parameters -= flex.double(self.components['decay'].n_params, maxB)
      self.components['decay'].calculate_scales_and_derivatives()
      logger.info('The "decay" model component has been rescaled, so that the\n'
        'maximum B-factor applied to any reflection is 0.0.')

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
        configdict['n_time_param']), col_names=['normalised_res_values',
        'norm_time_values'], parameter_esds=decay_setup['parameter_esds'])})
    if 'absorption' in configdict['corrections']:
      abs_setup = parameters_dict['absorption']
      self._components.update({'absorption' : SmoothScaleComponent3D(
        abs_setup['parameters'], shape=(configdict['n_x_param'],
        configdict['n_y_param'], configdict['n_time_param']),
        col_names=['normalised_x_abs_values', 'normalised_y_abs_values',
          'norm_time_values'],
        parameter_esds=abs_setup['parameter_esds'])})
    if 'modulation' in configdict['corrections']:
      mod_setup = parameters_dict['modulation']
      self._components.update({'modulation' : SmoothScaleComponent2D(
        mod_setup['parameters'], shape=(configdict['n_x_mod_param'],
        configdict['n_y_mod_param']), col_names=['normalised_x_det_values',
        'normalised_y_det_values'], parameter_esds=mod_setup['parameter_esds'])})

  @property
  def consecutive_scaling_order(self):
    return [['decay'], ['absorption'], ['modulation']]

  def configure_reflection_table(self, reflection_table, experiment, params):
    refl_table = reflection_table
    refl_table['norm_time_values'] = (refl_table['xyzobs.px.value'].parts()[2]
        * self.configdict['time_norm_fac'])
    if 'decay' in self.components:
      d0_sel = refl_table['d'] == 0.0
      refl_table['d'].set_selected(d0_sel, 1.0)  #set for now, then set back to zero later
      refl_table['normalised_res_values'] = (((1.0 / (refl_table['d']**2))
        - self.configdict['resmin']) / self.configdict['res_bin_width'])
      refl_table['normalised_res_values'].set_selected(d0_sel, 0.0001)
      refl_table['d'].set_selected(d0_sel, 0.0)
    if 'absorption' in self.components:
      refl_table['normalised_x_abs_values'] = ((refl_table['xyzobs.px.value'].parts()[0]
        - self.configdict['xmin']) / self.configdict['x_bin_width'])
      refl_table['normalised_y_abs_values'] = ((refl_table['xyzobs.px.value'].parts()[1]
        - self.configdict['ymin']) / self.configdict['y_bin_width'])
    if 'modulation' in self.components:
      refl_table['normalised_x_det_values'] = ((refl_table['xyzobs.px.value'].parts()[0]
        - self.configdict['xmin']) / self.configdict['x_det_bin_width'])
      refl_table['normalised_y_det_values'] = ((refl_table['xyzobs.px.value'].parts()[1]
        - self.configdict['ymin']) / self.configdict['y_det_bin_width'])
    return refl_table

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

  @property
  def consecutive_scaling_order(self):
    return [['scale', 'decay']]

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
