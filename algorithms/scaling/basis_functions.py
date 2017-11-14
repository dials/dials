'''
Classes that take in a data manager object and minimisation parameters and
return the scale factors and derivatives of the scale factors w.r.t.
the parameters
'''
from cctbx.array_family import flex
import numpy as np

class basis_function(object):
  '''Superclass for basis function that takes in a data manager object and
  minimisation parameters and calcuates the scale factors and derivatives'''
  def __init__(self, data_manager_object, apm):
    self.data_manager = data_manager_object
    self.apm = apm

  def calculate_scale_factors(self):
    factors_to_multiply = []
    for i, active_param in enumerate(self.apm.active_parameterisation):
      SF_object = self.data_manager.g_parameterisation[active_param]
      SF_object.set_scale_factors(self.apm.x[self.apm.cumulative_active_params[i]:
                                                  self.apm.cumulative_active_params[i+1]])
      factors_to_multiply.append(self.data_manager.g_parameterisation[active_param].calculate_smooth_scales())
    return flex.double(np.prod(np.array(factors_to_multiply), axis=0)) * self.apm.constant_g_values
    
  def calculate_derivatives(self):
    if len(self.data_manager.g_parameterisation) == 1:
      #only one active parameter, so don't need to chain rule any derivatives
      active_param = self.data_manager.g_parameterisation.values()[0]
      return active_param.calculate_smooth_derivatives()
    else:
      derivatives = flex.double([])
      for i, active_param in enumerate(self.apm.active_parameterisation):
        derivs = self.data_manager.g_parameterisation[active_param].calculate_smooth_derivatives()
        scale_multipliers = []
        for param in self.data_manager.g_parameterisation.keys():
          if param != active_param:
            scale_multipliers.append(self.data_manager.g_parameterisation[param].get_scales_of_reflections())
        scale_multiplier = flex.double(np.prod(np.array(scale_multipliers), axis=0))
        tiling_factor = self.apm.cumulative_active_params[i+1] - self.apm.cumulative_active_params[i]
        derivatives.extend(derivs * flex.double(np.tile(scale_multiplier, tiling_factor)))
      return derivatives

  def return_basis(self):
    '''Return the calculated scale factors and derivatives to be used
    in minimisation'''
    return self.calculate_scale_factors(), self.calculate_derivatives()

class KB_basis_function(basis_function):
  def calculate_scale_factors(self):
    print "Current scale and B-factor are %s" % (list(self.apm.x))
    scale_term = self.apm.x[0:1]
    decay_term = self.apm.x[1:2]
    self.data_manager.g_scale.set_scale_factors(scale_term)
    self.data_manager.g_decay.set_scale_factors(decay_term)
    d = self.data_manager.g_decay.d_values
    return (flex.double([scale_term[0]]*len(d)) * 
            flex.double(np.exp(flex.double([decay_term[0]]*len(d))/(2.0*(d**2)))))

  def calculate_derivatives(self):
    d = self.data_manager.g_decay.d_values
    B = self.data_manager.g_decay.scale_factors
    decay_term = flex.double(np.exp(flex.double([B[0]]*len(d))/(2.0*(d**2))))
    scale_term = self.data_manager.g_scale.scale_factors
    derivatives = flex.double([])
    derivatives.extend(decay_term)
    derivatives.extend(flex.double([scale_term[0]]*len(d)) * decay_term / (2*(d**2)))
    return derivatives

class xds_basis_function_log(basis_function):
  '''Subclass of basis_function for xds parameterisation'''
  def calculate_scale_factors(self):
    SF_object = self.data_manager.g_parameterisation[self.data_manager.active_parameterisation]
    SF_object.set_scale_factors(self.apm.x)
    gxvalues = SF_object.calculate_smooth_scales()
    scale_factors = gxvalues + self.data_manager.constant_g_values
    return flex.double(np.exp(scale_factors))

  def calculate_derivatives(self):
    '''Derivatives are fixed by the parameterisation'''
    derivatives = self.data_manager.g_parameterisation[
      self.data_manager.active_parameterisation].calculate_smooth_derivatives()
    return derivatives
