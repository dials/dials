'''
Collection of factories for creating the scaling models.
To add a new scaling model, one must define a new extension
in dials.extensions.scaling_model_ext, create a new factory
in this file and a new model is dials.algorithms.scaling.model.
'''
from dials.array_family import flex
from libtbx.phil import parse
import dials.algorithms.scaling.Model as Model
import pkg_resources

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

class AimlessSMFactory(object):
  '''
  Factory for creating an aimless-like scaling model.
  '''
  @classmethod
  def create(cls, params, experiments, reflections):
    '''
    create the scaling model defined by the params.
    '''
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

    return Model.AimlessScalingModel(s_norm_fac, scale_parameters, d_norm_fac,
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
  def create(cls, _, __, ____):
    '''create the simple KB scaling model.'''
    K_param = flex.double([1.0])
    B_param = flex.double([0.0])
    return Model.KBScalingModel(K_param, B_param)
