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

class xds_basis_function(basis_function):
  '''Subclass of basis_function for xds parameterisation'''
  def calculate_scale_factors(self):
    SF_object = self.data_manager.g_parameterisation[self.data_manager.active_parameterisation]['object']
    SF_object.set_scale_factors(self.parameters)
    gxvalues = SF_object.calculate_smooth_scales()
    scale_factors = gxvalues * self.data_manager.constant_g_values
    return scale_factors

  def calculate_derivatives(self):
    '''Derivatives are fixed by the parameterisation'''
    derivatives = self.data_manager.g_parameterisation[
      self.data_manager.active_parameterisation]['object'].calculate_smooth_derivatives()
    return derivatives

class aimless_basis_function(basis_function):
  '''Subclass of basis_function for aimless parameterisation'''
  def calculate_scale_factors(self):
    ngscale = self.data_manager.n_g_scale_params
    ngdecay = self.data_manager.n_g_decay_params
    self.data_manager.g_scale.set_scale_factors(self.parameters[0:ngscale])
    self.data_manager.g_decay.set_scale_factors(self.parameters[ngscale:
                                                                ngscale+ngdecay])
    scale = self.data_manager.g_scale.calculate_smooth_scales()
    B = self.data_manager.g_decay.calculate_smooth_scales()
    B_factor = flex.double(np.exp(B/(2.0 * (
      self.data_manager.reflections_for_scaling['d'] **2))))
    self.data_manager.reflections_for_scaling['B_factor'] = B_factor
    self.data_manager.reflections_for_scaling['scale_factor'] = scale
    return scale * B_factor

  def calculate_derivatives(self):
    scale_derivatives = self.data_manager.g_scale.calculate_smooth_derivatives()
    B_derivatives = self.data_manager.g_decay.calculate_smooth_derivatives()
    d = self.data_manager.reflections_for_scaling['d']
    B = self.data_manager.g_decay.get_scales_of_reflections()
    B_factor = flex.double(np.exp(B/(2.0 * (d**2))))
    dTdB = (flex.double(np.exp(B/(2.0 * (d**2)))) *
            self.data_manager.g_scale.get_scales_of_reflections() / (2.0 * (d**2)))
    B_factor_derivatives = flex.double(np.tile(dTdB, self.data_manager.n_g_decay_params))
    derivatives = flex.double([])
    derivatives.extend(scale_derivatives * flex.double(
      np.tile(B_factor, self.data_manager.n_g_scale_params)))
    derivatives.extend(B_factor_derivatives * B_derivatives)
    return derivatives

class xds_basis_function_log(basis_function):
  '''Subclass of basis_function for xds parameterisation'''
  def calculate_scale_factors(self):
    gxvalues = flex.double([self.parameters[i] for i in
           self.data_manager.reflections_for_scaling[self.data_manager.active_bin_index]])
    scale_factors = gxvalues + self.data_manager.constant_g_values
    return flex.double(np.exp(scale_factors))

  def calculate_derivatives(self):
    '''Derivatives are fixed by the parameterisation'''
    derivatives = self.data_manager.g_parameterisation[
      self.data_manager.active_parameterisation]['derivatives']
    return derivatives
