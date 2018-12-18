"""
Definitions of scaling models - collections of scale components with appropriate
methods to define how these are composed into one model.
"""
from __future__ import print_function

import abc
import logging
from collections import OrderedDict
import numpy as np
from dials.array_family import flex
from dials.algorithms.scaling.model.components.scale_components import \
  SingleScaleFactor, SingleBScaleFactor, SHScaleComponent
from dials.algorithms.scaling.model.components.smooth_scale_components import \
  SmoothScaleComponent1D, SmoothBScaleComponent1D, SmoothScaleComponent2D,\
  SmoothScaleComponent3D
from dials.algorithms.scaling.scaling_utilities import sph_harm_table

logger = logging.getLogger('dials')

class ScalingModelBase(object):
  """Base class for scaling models."""

  id_ = None

  __metaclass__ = abc.ABCMeta

  def __init__(self, configdict, is_scaled=False):
    self._components = OrderedDict()
    self._configdict = configdict
    self._is_scaled = is_scaled
    self._error_model = None

  @property
  def is_scaled(self):
    """Indictor as to whether this model has previously been refined."""
    return self._is_scaled

  def limit_image_range(self, new_image_range):
    """Modify the model if necessary due to excluding batches."""
    pass

  def set_scaling_model_as_scaled(self):
    """Indicate a scaling process has been performed on the data."""
    self._is_scaled = True

  def set_scaling_model_as_unscaled(self):
    """Indicate that no scaled data is associated with this model."""
    self._is_scaled = False

  def configure_reflection_table(self, reflection_table, experiment, params):
    """Perform calculations necessary to update the reflection table."""

  def set_valid_image_range(self, image_range):
    """Track the batch range for which the model corresponds to."""
    self._configdict['valid_image_range'] = image_range

  def normalise_components(self):
    """Optionally define a normalisation of the parameters after scaling."""

  @property
  def error_model(self):
    """The error model associated with the scaling model."""
    return self._error_model

  @property
  def configdict(self):
    """Dictionary of configuration parameters."""
    return self._configdict

  @property
  def components(self):
    """The components of the model, a dictionary."""
    return self._components

  @abc.abstractproperty
  def consecutive_refinement_order(self):
    """Return a nested list of correction names, to indicate the order
    to perform scaling in consecutive scaling mode if concurrent=0.
    e.g. [['scale', 'decay'], ['absorption']] would cause the first cycle to
    refine scale and decay, and then absorption in a subsequent cycle."""

  def to_dict(self):
    """Format data to dictionary for output."""
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
  @abc.abstractmethod
  def from_dict(cls, obj):
    """Create a scaling model object from a dictionary."""

  def set_error_model(self, error_model):
    """Associate an error model with the dataset."""
    self._error_model = error_model
    self._configdict.update({
      'error_model_type' : self.error_model.__class__.__name__,
      'error_model_parameters' : list(error_model.refined_parameters)})

  def show(self):
    """Print a representation of the scaling model."""
    print("Warning: Use of the .show() method is deprecated. Use print(object) instead.")
    print(str(self))

  def __str__(self):
    msg = ["Scaling model:"]
    msg.append("  type : " + str(self.id_))
    for name, component in self.components.iteritems():
      msg.append("  " + str(name).capitalize() + " component:")
      if component.parameter_esds:
        msg.append("    parameters (sigma)")
        for p, e in zip(component.parameters, component.parameter_esds):
          if p < 0.0:
            msg.append("    %.4f   (%.4f)" % (p, e))
          else:
            msg.append("     %.4f   (%.4f)" % (p, e))
      else:
        msg.append("    parameters")
        for p in component.parameters:
          if p < 0.0:
            msg.append("    %.4f" % p)
          else:
            msg.append("     %.4f" % p)
    msg.append("")
    return "\n".join(msg)

class PhysicalScalingModel(ScalingModelBase):
  """A scaling model for a physical parameterisation."""

  id_ = 'physical'

  def __init__(self, parameters_dict, configdict, is_scaled=False):
    super(PhysicalScalingModel, self).__init__(configdict, is_scaled)
    if 'scale' in configdict['corrections']:
      scale_setup = parameters_dict['scale']
      self._components.update({'scale' : SmoothScaleComponent1D(
        scale_setup['parameters'], scale_setup['parameter_esds'])})
    if 'decay' in configdict['corrections']:
      decay_setup = parameters_dict['decay']
      self._components.update({'decay' : SmoothBScaleComponent1D(
        decay_setup['parameters'], decay_setup['parameter_esds'])})
    if 'absorption' in configdict['corrections']:
      absorption_setup = parameters_dict['absorption']
      self._components.update({'absorption' : SHScaleComponent(
        absorption_setup['parameters'], absorption_setup['parameter_esds'])})

  @property
  def consecutive_refinement_order(self):
    return [['scale', 'decay'], ['absorption']]

  def configure_reflection_table(self, reflection_table, experiment, params):

    if 'scale' in self.components:
      norm = reflection_table['xyzobs.px.value'].parts()[2] * self._configdict['s_norm_fac']
      self.components['scale'].data = {'x': norm}
    if 'decay' in self.components:
      norm = reflection_table['xyzobs.px.value'].parts()[2] * self._configdict['d_norm_fac']
      self.components['decay'].parameter_restraints = flex.double(
        self.components['decay'].parameters.size(),
        params.parameterisation.decay_restraint)
      self.components['decay'].data = {'x' : norm, 'd' : reflection_table['d']}
    if 'absorption' in self.components:
      lmax = self._configdict['lmax']
      #here just pass in good reflections
      self.components['absorption'].data['sph_harm_table'] = sph_harm_table(
        reflection_table, experiment, lmax)
      surface_weight = self._configdict['abs_surface_weight']
      parameter_restraints = flex.double([])
      for i in range(1, lmax+1):
        parameter_restraints.extend(flex.double([1.0] * ((2*i)+1)))
      parameter_restraints *= surface_weight
      self.components['absorption'].parameter_restraints = parameter_restraints

  def limit_image_range(self, new_image_range):
    """Change the model to be suitable for a reduced batch range"""
    conf = self.configdict
    current_image_range = conf['valid_image_range']
    current_osc_range = conf['valid_osc_range']
    # calculate one osc as don't have access to scan object here
    one_osc = (conf['valid_osc_range'][1] - conf['valid_osc_range'][0])/(
        (conf['valid_image_range'][1] - (conf['valid_image_range'][0] - 1)))
    new_osc_range = ((new_image_range[0] - current_image_range[0]) * one_osc,
        (new_image_range[1] - current_image_range[0] + 1) * one_osc)
    if 'scale' in self.components:
      n_param, s_norm_fac, scale_rot_int = initialise_smooth_input(
        new_osc_range, one_osc, conf['scale_rot_interval'])
      n_old_params = len(self.components['scale'].parameters)
      conf['scale_rot_interval'] = scale_rot_int
      conf['s_norm_fac'] = s_norm_fac
      offset = calculate_new_offset(current_image_range[0], new_image_range[0],
        s_norm_fac, n_old_params, n_param)
      new_params = self.components['scale'].parameters[offset : offset + n_param]
      self.components['scale'].set_new_parameters(new_params)
    if 'decay' in self.components:
      n_param, d_norm_fac, decay_rot_int = initialise_smooth_input(
        new_osc_range, one_osc, conf['decay_rot_interval'])
      n_old_params = len(self.components['decay'].parameters)
      conf['decay_rot_interval'] = decay_rot_int
      conf['d_norm_fac'] = d_norm_fac
      offset = calculate_new_offset(current_image_range[0], new_image_range[0],
        d_norm_fac, n_old_params, n_param)
      new_params = self.components['decay'].parameters[offset : offset + n_param]
      self.components['decay'].set_new_parameters(new_params)

    new_osc_range_0 = current_osc_range[0] + (
      (new_image_range[0] - current_image_range[0]) * one_osc)
    new_osc_range_1 = current_osc_range[1] + (
      (new_image_range[1] - current_image_range[1]) * one_osc)
    self._configdict['valid_osc_range'] = (new_osc_range_0, new_osc_range_1)
    self.set_valid_image_range(new_image_range)

  def normalise_components(self):
    if 'scale' in self.components:
      # Do an invariant rescale of the scale at t=0 to one.'''
      joined_norm_vals = flex.double([])
      joined_inv_scales = flex.double([])
      for i in range(len(self.components['scale'].normalised_values)):
        scales, _ = self.components['scale'].calculate_scales_and_derivatives(block_id=i)
        norm_vals = self.components['scale'].normalised_values[i]
        joined_norm_vals.extend(norm_vals)
        joined_inv_scales.extend(scales)
      sel = (joined_norm_vals == min(joined_norm_vals))
      initial_scale = joined_inv_scales.select(sel)[0]
      self.components['scale'].parameters /= initial_scale
      logger.info('\nThe "scale" model component has been rescaled, so that the\n'
        'initial scale is 1.0.')
    if 'decay' in self.components:
      # Do an invariant rescale of the max B to zero.'''
      joined_d_vals = flex.double([])
      joined_inv_scales = flex.double([])
      for i in range(len(self.components['decay'].d_values)):
        scales, _ = self.components['decay'].calculate_scales_and_derivatives(block_id=i)
        d_vals = self.components['decay'].d_values[i]
        joined_d_vals.extend(d_vals)
        joined_inv_scales.extend(scales)
      maxB = flex.max(flex.double(np.log(joined_inv_scales))
                  * 2.0 * (joined_d_vals**2))
      self.components['decay'].parameters -= flex.double(
        self.components['decay'].n_params, maxB)
      logger.info('The "decay" model component has been rescaled, so that the\n'
        'maximum B-factor applied to any reflection is 0.0.')

  @classmethod
  def from_dict(cls, obj):
    """Create a scaling model object from a dictionary."""
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
  """A scaling model for an array-based parameterisation."""

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

  @property
  def consecutive_refinement_order(self):
    return [['decay'], ['absorption'], ['modulation']]

  def configure_reflection_table(self, reflection_table, _, __):
    xyz = reflection_table['xyzobs.px.value'].parts()
    norm_time = (xyz[2] * self.configdict['time_norm_fac'])
    if 'decay' in self.components:
      d = reflection_table['d']
      norm_res = (((1.0 / (d**2)) - self.configdict['resmin'])
        / self.configdict['res_bin_width'])
      self.components['decay'].data = {'x' : norm_res, 'y': norm_time}
    if 'absorption' in self.components:
      norm_x_abs = ((xyz[0] - self.configdict['xmin']) /
        self.configdict['x_bin_width'])
      norm_y_abs = ((xyz[1] - self.configdict['ymin']) /
        self.configdict['y_bin_width'])
      self.components['absorption'].data = {'x' : norm_x_abs,
        'y' : norm_y_abs, 'z': norm_time}
    if 'modulation' in self.components:
      norm_x_det = ((xyz[0] - self.configdict['xmin']) /
        self.configdict['x_det_bin_width'])
      norm_y_det = ((xyz[1] - self.configdict['ymin']) /
        self.configdict['y_det_bin_width'])
      self.components['modulation'].data = {'x' : norm_x_det, 'y': norm_y_det}

  ##FIXME update to not use reflection table
  def limit_image_range(self, new_image_range):
    """Change the model to be suitable for a reduced batch range"""
    conf = self.configdict
    current_image_range = conf['valid_image_range']
    current_osc_range = conf['valid_osc_range']
    # calculate one osc as don't have access to scan object here
    one_osc = (conf['valid_osc_range'][1] - conf['valid_osc_range'][0])/(
        (conf['valid_image_range'][1] - (conf['valid_image_range'][0] - 1)))
    new_osc_range = ((new_image_range[0] - current_image_range[0]) * one_osc,
        (new_image_range[1] - current_image_range[0] + 1) * one_osc)

    n_param, time_norm_fac, time_rot_int = initialise_smooth_input(
        new_osc_range, one_osc, conf['time_rot_interval'])
    n_old_time_params = int(len(self.components['decay'].parameters)/
      self.components['decay'].n_x_params)
    offset = calculate_new_offset(current_image_range[0], new_image_range[0],
        time_norm_fac, n_old_time_params, n_param)

    if 'absorption' in self.components:
      params = self.components['absorption'].parameters
      n_x_params = self.components['absorption'].n_x_params
      n_y_params = self.components['absorption'].n_y_params
      #can't do simple slice as 3-dim array
      time_offset = offset * n_x_params * n_y_params
      new_params = params[time_offset : time_offset + \
        (n_param * n_x_params * n_y_params)]
      self.components['absorption'].set_new_parameters(new_params,
        shape=(n_x_params, n_y_params, n_param))
    if 'decay' in self.components:
      params = self.components['decay'].parameters
      n_decay_params = self.components['decay'].n_x_params
      #can't do simple slice as 2-dim array
      decay_offset = offset * n_decay_params
      new_params = params[decay_offset : decay_offset + \
        (n_param * n_decay_params)]
      self.components['decay'].set_new_parameters(new_params,
        shape=(n_decay_params, n_param))
    self._configdict['n_time_param'] = n_param
    self._configdict['time_norm_fac'] = time_norm_fac
    self._configdict['time_rot_interval'] = time_rot_int
    new_osc_range_0 = current_osc_range[0] + (
      (new_image_range[0] - current_image_range[0]) * one_osc)
    new_osc_range_1 = current_osc_range[1] + (
      (new_image_range[1] - current_image_range[1]) * one_osc)
    self._configdict['valid_osc_range'] = (new_osc_range_0, new_osc_range_1)
    self.set_valid_image_range(new_image_range)

  @classmethod
  def from_dict(cls, obj):
    """Create a scaling model object from a dictionary."""
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
  """A scaling model for a KB parameterisation."""

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

  def configure_reflection_table(self, reflection_table, experiment, params):
    if 'scale' in self.components:
      self.components['scale'].data = {'id': reflection_table['id']}
    if 'decay' in self.components:
      self.components['decay'].data = {'d': reflection_table['d']}

  @property
  def consecutive_refinement_order(self):
    return [['scale', 'decay']]

  @classmethod
  def from_dict(cls, obj):
    """Create a scaling model object from a dictionary."""
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

def map_old_to_new_range(old_range, new_range):
  """Calculate the offset and number of params needed for the new dataset"""
  offset = 0
  n_param = range_to_n_param(new_range)
  n_old_param = range_to_n_param(old_range)
  if n_old_param == n_param: #would only work for using function to reduce range
    return 0, n_param
  if new_range[0] > 0.0:
    #might need to shift
    if n_param > 3 and n_old_param > 3:
      offset = int((new_range[0] - 0.5) //1) + 1
    elif n_param < 4 and n_old_param > 3:
      # want to use parameter nearest to new_range[0]
      if new_range[0] % 1.0 < 0.5:
        offset = int((new_range[0] - 0.5) //1) + 2
      else:
        offset = int((new_range[0] - 0.5) //1) + 1
    else:
      if new_range[0] % 1.0 < 0.5:
        offset = int(new_range[0]//1)
      else:
        offset = int(new_range[0] //1) + 1
  return offset, n_param

def calculate_new_offset(old_batch_0, new_batch_0, new_norm_fac, n_old_param,
  n_new_param):
  if n_old_param == 2:
    return 0 #cant have less than two params
  batch_difference = (new_batch_0 - old_batch_0) * new_norm_fac
  n_to_shift = int(batch_difference // 1)
  if batch_difference % 1 > 0.5:
    n_to_shift += 1
  return min(n_old_param - n_new_param, n_to_shift) #cant shift by more
  #than difference between old and new

def range_to_n_param(range_):
  """Calculate the number of smoother parameters from a range"""
  if (range_[1] - range_[0]) < 1.0:
    n_param = 2
  elif (range_[1] - range_[0]) < 2.0:
    n_param = 3
  else:
    n_param = max(int(range_[1] - range_[0])+1, 3) + 2
  return n_param

def initialise_smooth_input(osc_range, one_osc_width, interval):
  """Calculate the number of parameters and norm_fac/rot_int."""
  interval += 0.00001
  if (osc_range[1] - osc_range[0]) < (2.0 * interval):
    if (osc_range[1] - osc_range[0]) <= interval:
      rot_int = osc_range[1] - osc_range[0]
      n_param = 2
    else:
      rot_int = ((osc_range[1] - osc_range[0])/2.0)
      n_param = 3
  else:
    n_bins = max(int((osc_range[1] - osc_range[0])/ interval)+1, 3)
    rot_int = (osc_range[1] - osc_range[0])/float(n_bins)
    n_param = n_bins + 2
  norm_fac = 0.999 * one_osc_width / rot_int #to make sure normalise values
  #fall within range of smoother.
  return n_param, norm_fac, rot_int
