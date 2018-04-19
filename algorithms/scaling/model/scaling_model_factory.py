'''
Collection of factories for creating the scaling models.
To add a new scaling model, one must define a new extension
in dials.extensions.scaling_model_ext, create a new factory
in this file and create a new model in dials.algorithms.scaling.model.
'''
from collections import OrderedDict
from dials.array_family import flex
import dials.algorithms.scaling.model.model as Model

class KBSMFactory(object):
  '''
  Factory for creating a KB scaling model.
  '''
  @classmethod
  def create(cls, params, __, ____):
    '''create the simple KB scaling model.'''
    corrections = []
    if params.parameterisation.decay_term:
      corrections.append('decay')
    if params.parameterisation.scale_term:
      corrections.append('scale')
    configdict = OrderedDict({'corrections': corrections})

    K_param = flex.double([1.0])
    B_param = flex.double([0.0])

    parameters_dict = {
      'scale': {'parameters' : K_param, 'parameter_esds' : None},
      'decay': {'parameters' : B_param, 'parameter_esds' : None}}

    return Model.KBScalingModel(parameters_dict, configdict)

class PhysicalSMFactory(object):
  """
  Factory for creating a physical scaling model.
  """

  @classmethod
  def create(cls, params, experiments, reflections):
    """Create the scaling model defined by the params."""

    configdict = OrderedDict({'corrections':[]})
    parameters_dict = {}

    osc_range = osc_range_check_for_user_excluded(experiments, reflections)
    one_osc_width = experiments.scan.get_oscillation()[1]

    if params.parameterisation.scale_term:
      configdict['corrections'].append('scale')
      n_scale_param, s_norm_fac, scale_rot_int = initialise_smooth_input(
        osc_range, one_osc_width, params.parameterisation.scale_interval)
      configdict.update({'s_norm_fac' : s_norm_fac,
        'scale_rot_interval' : scale_rot_int})
      parameters_dict['scale'] = {'parameters' : flex.double(n_scale_param, 1.0),
        'parameter_esds' : None}

    if params.parameterisation.decay_term:
      configdict['corrections'].append('decay')
      n_decay_param, d_norm_fac, decay_rot_int = initialise_smooth_input(
        osc_range, one_osc_width, params.parameterisation.decay_interval)
      configdict.update({'d_norm_fac' : d_norm_fac,
        'decay_rot_interval' : decay_rot_int})
      parameters_dict['decay'] = {'parameters' : flex.double(n_decay_param, 0.0),
        'parameter_esds' : None}

    if params.parameterisation.absorption_term:
      configdict['corrections'].append('absorption')
      lmax = params.parameterisation.lmax
      n_abs_param = (2*lmax) + (lmax**2)  #arithmetic sum formula (a1=3, d=2)
      configdict.update({'lmax' : lmax})
      parameters_dict['absorption'] = {'parameters' : flex.double(
        n_abs_param, 0.0), 'parameter_esds' : None}

    return Model.PhysicalScalingModel(parameters_dict, configdict)


def initialise_smooth_input(osc_range, one_osc_width, interval):
  """Calculate the number of parameters and norm_fac/rot_int."""
  interval += 0.001
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
  norm_fac = 0.9999 * one_osc_width / rot_int #to make sure normalise values
  #fall within range of smoother.
  return n_param, norm_fac, rot_int

def calc_n_param_from_bins(value_min, value_max, n_bins):
  """Return the correct number of bins for initialising the gaussian
  smoothers."""
  assert n_bins > 0
  assert isinstance(n_bins, int)
  bin_width = (value_max - value_min) / n_bins
  if n_bins == 1:
    n_param = 2
  elif n_bins == 2:
    n_param = 3
  else:
    n_param = n_bins + 2
  return n_param, bin_width

def osc_range_check_for_user_excluded(experiments, reflections):
  """Determine the oscillation range, allowing for user excluded range."""
  osc_range = experiments.scan.get_oscillation_range()
  one_osc_width = experiments.scan.get_oscillation()[1]
  user_excluded = reflections.get_flags(
    reflections.flags.user_excluded_in_scaling)
  if user_excluded.count(True) > 0:
    reflections_for_scaling = reflections.select(~user_excluded)
    reflections_for_scaling = reflections_for_scaling.select(
      reflections_for_scaling.get_flags(reflections_for_scaling.flags.integrated))
    max_osc = (max(reflections_for_scaling['xyzobs.px.value'].parts()[2]
      * one_osc_width) + experiments.scan.get_oscillation()[0])
    min_osc = (min(reflections_for_scaling['xyzobs.px.value'].parts()[2]
      * one_osc_width) + experiments.scan.get_oscillation()[0])
    if max_osc < osc_range[1] - 0.8: #some end frames excluded
      min_osc = osc_range[0]
      osc_range = (min_osc, max_osc + 0.001)
    elif min_osc > osc_range[0] + 0.8: #some beginning frames excluded
      max_osc = osc_range[1]
      osc_range = (min_osc, max_osc)
  return osc_range

class ArraySMFactory(object):
  """
  Factory for creating an array-based scaling model.
  """

  @classmethod
  def create(cls, params, experiments, reflections):
    '''create an array-based scaling model.'''
    reflections = reflections.select(reflections['d'] > 0.0)

    # First initialise things common to more than one correction.
    osc_range = osc_range_check_for_user_excluded(experiments, reflections)
    one_osc_width = experiments.scan.get_oscillation()[1]
    n_time_param, time_norm_fac, time_rot_int = initialise_smooth_input(
      osc_range, one_osc_width, params.parameterisation.decay_interval)
    (xvalues, yvalues, _) = reflections['xyzobs.px.value'].parts()
    (xmax, xmin) = (max(xvalues) + 0.001, min(xvalues) - 0.001)
    (ymax, ymin) = (max(yvalues) + 0.001, min(yvalues) - 0.001)

    corrections = []
    if params.parameterisation.decay_term:
      corrections.append('decay')
      resmax = (1.0 / (min(reflections['d'])**2)) + 0.001
      resmin = (1.0 / (max(reflections['d'])**2)) - 0.001
      n_res_bins = params.parameterisation.n_resolution_bins
      n_res_param, res_bin_width = calc_n_param_from_bins(resmin, resmax,
        n_res_bins)
      dec_params = flex.double((n_time_param * n_res_param), 1.0)
    if params.parameterisation.absorption_term:
      corrections.append('absorption')
      nxbins = nybins = params.parameterisation.n_absorption_bins
      n_x_param, x_bin_width = calc_n_param_from_bins(xmin, xmax, nxbins)
      n_y_param, y_bin_width = calc_n_param_from_bins(ymin, ymax, nybins)
      abs_params = flex.double((n_x_param * n_y_param * n_time_param), 1.0)
    if params.parameterisation.modulation_term:
      corrections.append('modulation')
      nx_det_bins = ny_det_bins = params.parameterisation.n_modulation_bins
      n_x_mod_param, x_det_bin_width = calc_n_param_from_bins(
        xmin, xmax, nx_det_bins)
      n_y_mod_param, y_det_bin_width = calc_n_param_from_bins(
        ymin, ymax, ny_det_bins)
      mod_params = flex.double((n_x_mod_param * n_y_mod_param), 1.0)

    configdict = OrderedDict({'corrections': corrections})
    parameters_dict = {}

    if 'decay' in configdict['corrections']:
      configdict.update({'n_res_param': n_res_param, 'n_time_param': n_time_param,
        'resmin' : resmin, 'res_bin_width' : res_bin_width,
        'time_norm_fac' : time_norm_fac, 'time_rot_interval' : time_rot_int})
      parameters_dict['decay'] = {'parameters' : dec_params,
        'parameter_esds' : None}

    if 'absorption' in configdict['corrections']:
      configdict.update({'n_x_param' : n_x_param, 'n_y_param' : n_y_param,
        'xmin' : xmin, 'ymin' : ymin, 'x_bin_width' : x_bin_width,
        'y_bin_width' : y_bin_width, 'n_time_param': n_time_param,
        'time_norm_fac' : time_norm_fac, 'time_rot_interval' : time_rot_int
      })
      parameters_dict['absorption'] = {'parameters' : abs_params,
        'parameter_esds' : None}

    if 'modulation' in configdict['corrections']:
      configdict.update({
        'n_x_mod_param' : n_x_mod_param, 'n_y_mod_param' : n_y_mod_param,
        'xmin' : xmin, 'ymin' : ymin,
        'x_det_bin_width' : x_det_bin_width, 'y_det_bin_width' : y_det_bin_width
      })
      parameters_dict['modulation'] = {'parameters' : mod_params,
        'parameter_esds' : None}

    return Model.ArrayScalingModel(parameters_dict, configdict)
