'''
Collection of factories for creating the scaling models.
To add a new scaling model, one must define a new extension
in dials.extensions.scaling_model_ext, create a new factory
in this file and a new model is dials.algorithms.scaling.model.
'''
from dials.array_family import flex
import dials.algorithms.scaling.model.Model as Model
import pkg_resources
from collections import OrderedDict


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
      else:
        for entry_point in pkg_resources.iter_entry_points('dxtbx.scaling_model_ext'):
          if entry_point.name == params.scaling_model:
            #finds relevant extension in dials.extensions.scaling_model_ext
            factory = entry_point.load().factory()
            exp.scaling_model = factory.create(params, exp, refl)
    return experiments

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

    return Model.KBScalingModel(K_param, B_param, configdict)

class AimlessSMFactory(object):
  '''
  Factory for creating an aimless-like scaling model.
  '''
  @classmethod
  def create(cls, params, experiments, reflections):
    '''
    create the scaling model defined by the params.
    '''
    corrections = []
    #these names are used many time in the program as flags, create a better
    #way to strongly enforce inclusion of certain corrections?
    if params.parameterisation.scale_term:
      corrections.append('scale')
    if params.parameterisation.decay_term:
      corrections.append('decay')
    if params.parameterisation.absorption_term:
      corrections.append('absorption')

    configdict = OrderedDict({'corrections': corrections})

    scale_rot_int = params.parameterisation.scale_interval + 0.001
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

    decay_rot_int = params.parameterisation.decay_interval + 0.001
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

    if 'scale' in configdict['corrections']:
      configdict.update({'s_norm_fac' : s_norm_fac})
    if 'decay' in configdict['corrections']:
      configdict.update({'d_norm_fac' : d_norm_fac})
    if 'absorption' in configdict['corrections']:
      configdict.update({'lmax' : lmax})

    return Model.AimlessScalingModel(scale_parameters, decay_parameters,
      abs_parameters, configdict)

class XscaleSMFactory(object):
  '''
  Factory for creating an scale-like scaling model.
  '''
  @classmethod
  def create(cls, params, experiments, reflections):
    '''create an XScale scaling model.'''
    #assert 0, 'method not yet implemented'
    reflections = reflections.select(reflections['d'] > 0.0)

    #only create components that are specified in params?
    corrections = []
    if params.parameterisation.decay_term:
      corrections.append('decay')
    if params.parameterisation.absorption_term:
      corrections.append('absorption')
    if params.parameterisation.modulation_term:
      corrections.append('modulation')
    configdict = OrderedDict({'corrections': corrections})

    scale_rot_int = params.parameterisation.scale_interval + 0.001
    osc_range = experiments.scan.get_oscillation_range()
    nzbins = int((osc_range[1] - osc_range[0])/ scale_rot_int)
    if ((osc_range[1] - osc_range[0])/ scale_rot_int) % 1 < 0.33:
      scale_rot_int = (osc_range[1] - osc_range[0])/float(nzbins) + 0.001
      nzbins = int((osc_range[1] - osc_range[0])/ scale_rot_int)

    '''Bin the data into resolution and time 'z' bins'''
    (xvalues, yvalues, zvalues) = reflections['xyzobs.px.value'].parts()
    (xmax, xmin) = (max(xvalues) + 0.001, min(xvalues) - 0.001)
    (ymax, ymin) = (max(yvalues) + 0.001, min(yvalues) - 0.001)
    (zmax, zmin) = (max(zvalues) + 0.001, min(zvalues) - 0.001)
    resmax = (1.0 / (min(reflections['d'])**2)) + 0.001
    resmin = (1.0 / (max(reflections['d'])**2)) - 0.001

    n_res_bins = 20 #params.scaling_options.n_res_bins

    res_bin_width = (resmax - resmin) / n_res_bins
    time_bin_width = (zmax - zmin) / nzbins
    nres = ((1.0 / (reflections['d']**2)) - resmin) / res_bin_width
    n_res_param = int(max(nres)//1) - int(min(nres)//1) + 3 #for g_decay
    nt = (zvalues - zmin) / time_bin_width
    n_time_param = int(max(nt)//1) - int(min(nt)//1) + 3 #for g_decay and g_abs
    dec_params = flex.double([1.0] * n_time_param * n_res_param)

    if 'decay' in configdict['corrections']:
      configdict.update({'n_res_param': n_res_param, 'n_time_param': n_time_param,
        'resmin' : resmin, 'res_bin_width' : res_bin_width, 'zmin' : zmin,
        'time_bin_width' : time_bin_width})

    nxbins = nybins = 3.0
    x_bin_width = (xmax - xmin) / float(nxbins)
    y_bin_width = (ymax - ymin) / float(nybins)
    nax = ((xvalues - xmin) / x_bin_width)
    n_x_param = int(max(nax)//1) - int(min(nax)//1) + 3
    nay = ((yvalues - ymin) / y_bin_width)
    n_y_param = int(max(nay)//1) - int(min(nay)//1) + 3
    abs_params = flex.double([1.0] * n_x_param * n_y_param * n_time_param)

    if 'absorption' in configdict['corrections']:
      configdict.update({
        'n_x_param' : n_x_param, 'n_y_param' : n_y_param, 'xmin' : xmin,
        'ymin' : ymin, 'x_bin_width' : x_bin_width, 'y_bin_width' : y_bin_width,
        'n_time_param': n_time_param, 'zmin' : zmin, 'time_bin_width' : time_bin_width
      })

    n_detector_bins = 10#params.parameterisation.n_detector_bins
    nx_det_bins = ny_det_bins = n_detector_bins
    x_det_bin_width = (xmax - xmin) / float(nx_det_bins)
    y_det_bin_width = (ymax - ymin) / float(ny_det_bins)
    nxdet = ((xvalues - xmin) / x_det_bin_width)
    n_x_mod_param = int(max(nxdet)//1) - int(min(nxdet)//1) + 3
    nydet = ((yvalues - ymin) / y_det_bin_width)
    n_y_mod_param = int(max(nydet)//1) - int(min(nydet)//1) + 3
    mod_params = flex.double([1.0] * n_x_mod_param * n_y_mod_param)

    if 'modulation' in configdict['corrections']:
      configdict.update({
        'n_x_mod_param' : n_x_mod_param, 'n_y_mod_param' : n_y_mod_param,
        'xmin' : xmin, 'ymin' : ymin,
        'x_det_bin_width' : x_det_bin_width, 'y_det_bin_width' : y_det_bin_width
      })

    return Model.XscaleScalingModel(dec_params, abs_params, mod_params, configdict)
