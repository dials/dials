'''
Classes that take in a data manager object and minimisation parameters and
return the scale factors and derivatives of the scale factors w.r.t.
the parameters
'''
from __future__ import print_function
#from cctbx.array_family import flex
from dials.array_family import flex
import numpy as np
from dials_scaling_helpers_ext import row_multiply
from scitbx import sparse

class basis_function(object):
  '''Superclass for basis function that takes in a data manager object and
  minimisation parameters and calcuates the scale factors and derivatives'''
  def __init__(self, data_manager_object, apm):
    self.data_manager = data_manager_object
    self.apm = apm

  def update_scale_factors(self):
    '''update the parameters in each SF object from the apm parameter list.'''
    for i, active_param in enumerate(self.apm.active_parameterisation):
      SF_object = self.data_manager.g_parameterisation[active_param]
      SF_object.parameters = self.apm.x[self.apm.cumulative_active_params[i]:
                                            self.apm.cumulative_active_params[i+1]]
      SF_object.calculate_scales_and_derivatives()

  def calculate_scale_factors(self):
    '''calculate overall scale factor from reflections from individual components'''
    multiplied_scale_factors = flex.double([1.0] * self.data_manager.Ih_table.size)
    for active_param in self.apm.active_parameterisation:
      multiplied_scale_factors *= self.data_manager.g_parameterisation[
        active_param].inverse_scales
    if self.apm.constant_g_values:
      multiplied_scale_factors *= self.apm.constant_g_values
    return multiplied_scale_factors

  def calculate_derivatives(self):
    '''calculate derivatives matrix'''
    if len(self.data_manager.g_parameterisation) == 1:
      #only one active parameter, so don't need to chain rule any derivatives
      active_param = self.data_manager.g_parameterisation.values()[0]
      return active_param.derivatives
    else:
      n_refl = self.data_manager.Ih_table.size
      derivatives = sparse.matrix(n_refl, self.apm.cumulative_active_params[-1])
      for i, active_param in enumerate(self.apm.active_parameterisation):
        derivs = self.data_manager.g_parameterisation[active_param].derivatives
        scale_multipliers = flex.double([1.0] * n_refl)
        for param, SF_obj in self.data_manager.g_parameterisation.iteritems():
          if param != active_param:
            scale_multipliers *= SF_obj.inverse_scales
        next_deriv = row_multiply(derivs, scale_multipliers)
        derivatives.assign_block(next_deriv, 0, self.apm.cumulative_active_params[i])
      return derivatives

  def return_basis(self):
    '''Return the calculated scale factors and derivatives to be used
    in minimisation'''
    self.update_scale_factors()
    return self.calculate_scale_factors(), self.calculate_derivatives()

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
