'''
Classes that take in a data manager object and minimisation parameters and
return the scale factors and derivatives of the scale factors w.r.t.
the parameters
'''
from __future__ import print_function
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
      factors_to_multiply.append(self.data_manager.g_parameterisation[
        active_param].calculate_smooth_scales())
    if self.apm.constant_g_values:
      return (flex.double(np.prod(np.array(factors_to_multiply), axis=0))
              * self.apm.constant_g_values)
    else:
      return flex.double(np.prod(np.array(factors_to_multiply), axis=0))
    
  def calculate_derivatives(self):
    if len(self.data_manager.g_parameterisation) == 1:
      #only one active parameter, so don't need to chain rule any derivatives
      active_param = self.data_manager.g_parameterisation.values()[0]
      return active_param.calculate_smooth_derivatives()
    else:
      derivatives = flex.double([])
      for i, active_param in enumerate(self.apm.active_parameterisation):
        derivs = self.data_manager.g_parameterisation[
          active_param].calculate_smooth_derivatives()
        scale_multipliers = []
        for param in self.data_manager.g_parameterisation.keys():
          if param != active_param:
            scale_multipliers.append(self.data_manager.g_parameterisation[
              param].get_scales_of_reflections())
        scale_mult = flex.double(np.prod(np.array(scale_multipliers), axis=0))
        tile_factor = (self.apm.cumulative_active_params[i+1]
                       - self.apm.cumulative_active_params[i])
        derivatives.extend(derivs * flex.double(np.tile(scale_mult, tile_factor)))
      return derivatives

  def return_basis(self):
    '''Return the calculated scale factors and derivatives to be used
    in minimisation'''
    return self.calculate_scale_factors(), self.calculate_derivatives()

class KB_basis_function(basis_function):
  def __init__(self, data_manager_object, apm):
    basis_function.__init__(self, data_manager_object, apm)

  def calculate_scale_factors(self):
    factors_to_multiply = []
    msg = ('Current parameters: {0}= {1}').format(
      ''.join(i.lstrip('g_')+' ' for i in self.data_manager.g_parameterisation),
      ' '.join('%.6f' % i for i in list(self.apm.x)))
    print(msg)
    for i, active_param in enumerate(self.apm.active_parameterisation):
      SF_object = self.data_manager.g_parameterisation[active_param]
      SF_object.set_scale_factors(self.apm.x[self.apm.cumulative_active_params[i]:
                                             self.apm.cumulative_active_params[i+1]])
      factors_to_multiply.append(self.data_manager.g_parameterisation[
        active_param].get_scales_of_reflections())
    return flex.double(np.prod(np.array(factors_to_multiply), axis=0))

  def calculate_derivatives(self):
    if len(self.data_manager.g_parameterisation) == 1:
      #only one active parameter, so don't need to chain rule any derivatives
      active_param = self.data_manager.g_parameterisation.values()[0]
      return active_param.get_derivatives()
    else:
      derivatives = flex.double([])
      for i, active_param in enumerate(self.apm.active_parameterisation):
        derivs = self.data_manager.g_parameterisation[active_param].get_derivatives()
        scale_multipliers = []
        for param in self.data_manager.g_parameterisation.keys():
          if param != active_param:
            scale_multipliers.append(self.data_manager.g_parameterisation[
              param].get_scales_of_reflections())
        scale_mult = flex.double(np.prod(np.array(scale_multipliers), axis=0))
        tile_factor = (self.apm.cumulative_active_params[i+1]
                       - self.apm.cumulative_active_params[i])
        derivatives.extend(derivs * flex.double(np.tile(scale_mult, tile_factor)))
      return derivatives

class xds_basis_function_log(basis_function):
  '''Subclass of basis_function for xds parameterisation'''
  def calculate_scale_factors(self):
    SF_object = self.data_manager.g_parameterisation[
      self.data_manager.active_parameterisation]
    SF_object.set_scale_factors(self.apm.x)
    gxvalues = SF_object.calculate_smooth_scales()
    scale_factors = gxvalues + self.data_manager.constant_g_values
    return flex.double(np.exp(scale_factors))

  def calculate_derivatives(self):
    '''Derivatives are fixed by the parameterisation'''
    derivatives = self.data_manager.g_parameterisation[
      self.data_manager.active_parameterisation].calculate_smooth_derivatives()
    return derivatives
