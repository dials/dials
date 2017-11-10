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
  def __init__(self, data_manager_object, parameters=None):
    self.data_manager = data_manager_object
    self.parameters = parameters

  def calculate_scale_factors(self):
    '''Calculate scale factors method to be overloaded by subclass'''
    assert 0, 'parameterization-specific scale factors must be defined in a subclass'

  def calculate_derivatives(self):
    '''Calculate derivatives method to be overloaded by subclass'''
    assert 0, 'parameterization-specific derivatives must be defined in a subclass'

  def return_basis(self):
    '''Return the calculated scale factors and derivatives to be used
    in minimisation'''
    return self.calculate_scale_factors(), self.calculate_derivatives()

class KB_basis_function(basis_function):
  def calculate_scale_factors(self):
    print "Current scale and B-factor are %s" % (list(self.parameters))
    scale_term = self.parameters[0:1]
    decay_term = self.parameters[1:2]
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

class xds_basis_function(basis_function):
  '''Subclass of basis_function for xds parameterisation'''
  def calculate_scale_factors(self):
    SF_object = self.data_manager.g_parameterisation[self.data_manager.active_parameterisation]
    SF_object.set_scale_factors(self.parameters)
    gxvalues = SF_object.calculate_smooth_scales()
    scale_factors = gxvalues * self.data_manager.constant_g_values
    return scale_factors

  def calculate_derivatives(self):
    '''Derivatives are fixed by the parameterisation'''
    derivatives = self.data_manager.g_parameterisation[
      self.data_manager.active_parameterisation].calculate_smooth_derivatives() 
    #should this be multiplied by constant_g_values?
    return derivatives * flex.double(np.tile(self.data_manager.constant_g_values, 
                                             self.data_manager.n_active_params))

class aimless_basis_function(basis_function):
  '''Subclass of basis_function for aimless parameterisation'''
  def calculate_scale_factors(self):
    ngscale = self.data_manager.n_g_scale_params
    ngdecay = self.data_manager.n_g_decay_params
    #ngabs = self.data_manager.n_g_abs_params
    factors_to_multiply = []
    for i, active_param in enumerate(self.data_manager.active_parameterisation):
      self.data_manager.g_parameterisation[active_param].set_scale_factors(
        self.parameters[self.data_manager.cumulative_active_params[i]:self.data_manager.cumulative_active_params[i+1]])
      factors_to_multiply.append(self.data_manager.g_parameterisation[active_param].calculate_smooth_scales())
    return flex.double(np.prod(np.array(factors_to_multiply), axis=0)) * self.data_manager.constant_g_values
    
    '''if self.data_manager.scaling_options['scale_term']:
      self.data_manager.g_scale.set_scale_factors(self.parameters[0:ngscale])
      factors_to_multiply.append(self.data_manager.g_scale.calculate_smooth_scales())
    if self.data_manager.scaling_options['decay_term']:
      self.data_manager.g_decay.set_scale_factors(self.parameters[ngscale:
                                                                  ngscale+ngdecay])
      factors_to_multiply.append(self.data_manager.g_decay.calculate_smooth_scales())
    if self.data_manager.scaling_options['absorption_term']:
      self.data_manager.g_absorption.set_scale_factors(self.parameters[ngscale+ngdecay:])
      factors_to_multiply.append(self.data_manager.g_absorption.calculate_scales())'''

    #return flex.double(np.prod(np.array(factors_to_multiply), axis=0))
    #scale = self.data_manager.g_scale.calculate_smooth_scales()
    #B = self.data_manager.g_decay.calculate_smooth_scales()
    #S = self.data_manager.g_absorption.calculate_scales()
    #return scale * B * S

  def calculate_derivatives(self):
    derivatives = flex.double([])
    for i, active_param in enumerate(self.data_manager.active_parameterisation):
      derivs = self.data_manager.g_parameterisation[active_param].calculate_smooth_derivatives()
      scale_multipliers = []
      for param in ['g_scale', 'g_absorption', 'g_decay']:
        if param != active_param:
          scale_multipliers.append(self.data_manager.g_parameterisation[param].get_scales_of_reflections())
      scale_multiplier = flex.double(np.prod(np.array(scale_multipliers), axis=0))
      tiling_factor = self.data_manager.cumulative_active_params[i+1] - self.data_manager.cumulative_active_params[i]
      derivatives.extend(derivs * flex.double(np.tile(scale_multiplier, tiling_factor)))
    return derivatives

    '''scale_derivatives = self.data_manager.g_scale.calculate_smooth_derivatives()
    B_derivatives = self.data_manager.g_decay.calculate_smooth_derivatives()
    d = self.data_manager.g_decay.d_values#reflections_for_scaling['d']
    B = self.data_manager.g_decay.get_scales_of_reflections()
    S = self.data_manager.g_absorption.get_scales_of_reflections()
    scale_term = self.data_manager.g_scale.get_scales_of_reflections()
    #B_factor = flex.double(np.exp(B/(2.0 * (d**2))))
    dTdB = (B * scale_term * S / (2.0 * (d**2)))
    B_factor_derivatives = flex.double(np.tile(dTdB, self.data_manager.n_g_decay_params))
    derivatives = flex.double([])
    derivatives.extend(scale_derivatives * flex.double(
      np.tile((B * S), self.data_manager.n_g_scale_params)))
    derivatives.extend(B_factor_derivatives * B_derivatives)
    if self.data_manager.scaling_options['absorption_term']:
      abs_derivatives = self.data_manager.g_absorption.derivatives
      derivatives.extend(abs_derivatives * flex.double(
        np.tile((B * scale_term), self.data_manager.n_g_abs_params)))
    return derivatives'''

class xds_basis_function_log(basis_function):
  '''Subclass of basis_function for xds parameterisation'''
  def calculate_scale_factors(self):
    SF_object = self.data_manager.g_parameterisation[self.data_manager.active_parameterisation]
    SF_object.set_scale_factors(self.parameters)
    gxvalues = SF_object.calculate_smooth_scales()
    scale_factors = gxvalues + self.data_manager.constant_g_values
    return flex.double(np.exp(scale_factors))

  def calculate_derivatives(self):
    '''Derivatives are fixed by the parameterisation'''
    derivatives = self.data_manager.g_parameterisation[
      self.data_manager.active_parameterisation].calculate_smooth_derivatives()
    return derivatives
