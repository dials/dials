'''
Classes that take in a scaler and minimisation parameters and
return the scale factors and derivatives of the scale factors w.r.t.
the parameters
'''
from __future__ import print_function
from dials.array_family import flex
from dials_scaling_helpers_ext import row_multiply
from scitbx import sparse

class basis_function(object):
  '''Superclass for basis function that takes in a scaler and
  minimisation parameters and calcuates the scale factors and derivatives'''
  def __init__(self, scaler, apm):
    self.scaler = scaler
    self.apm = apm

  def update_scale_factors(self):
    '''update the parameters in each SF object from the apm parameter list.'''
    for i, active_param in enumerate(self.apm.active_parameterisation):
      SF_object = self.scaler.experiments.scaling_model.components[active_param]
      SF_object.parameters = self.apm.x[self.apm.n_cumul_params_list[i]:
        self.apm.n_cumul_params_list[i+1]]
      SF_object.calculate_scales_and_derivatives()

  def calculate_scale_factors(self):
    '''calculate overall scale factor from reflections from individual components'''
    multiplied_scale_factors = flex.double([1.0] * self.scaler.Ih_table.size)
    for active_param in self.apm.active_parameterisation:
      multiplied_scale_factors *= self.scaler.experiments.scaling_model.components[
        active_param].inverse_scales
    if self.apm.constant_g_values:
      multiplied_scale_factors *= self.apm.constant_g_values
    return multiplied_scale_factors

  def calculate_derivatives(self):
    '''calculate derivatives matrix'''
    if len(self.scaler.experiments.scaling_model.components) == 1:
      #only one active parameter, so don't need to chain rule any derivatives
      return self.scaler.experiments.scaling_model.components.values()[0].derivatives
    derivatives = sparse.matrix(self.scaler.Ih_table.size,
      self.apm.n_active_params)
    for i, active_param in enumerate(self.apm.active_parameterisation):
      derivs = self.scaler.experiments.scaling_model.components[active_param].derivatives
      scale_multipliers = flex.double([1.0] * self.scaler.Ih_table.size)
      for param, SF_obj in self.scaler.experiments.scaling_model.components.iteritems():
        if param != active_param:
          scale_multipliers *= SF_obj.inverse_scales
      next_deriv = row_multiply(derivs, scale_multipliers)
      derivatives.assign_block(next_deriv, 0, self.apm.n_cumul_params_list[i])
    return derivatives

  def return_basis(self):
    '''Return the calculated scale factors and derivatives to be used
    in minimisation'''
    self.update_scale_factors()
    return self.calculate_scale_factors(), self.calculate_derivatives()

'''class xds_basis_function_log(basis_function):
  'Subclass of basis_function for logarithmic parameterisation'
  def calculate_scale_factors(self):
    SF_object = self.scaler.g_parameterisation[self.scaler.active_parameterisation]
    SF_object.set_scale_factors(self.apm.x)
    gxvalues = SF_object.calculate_smooth_scales()
    scale_factors = gxvalues + self.scaler.constant_g_values
    return flex.double(np.exp(scale_factors))

  def calculate_derivatives(self):
    'Derivatives are fixed by the parameterisation'
    derivatives = self.scaler.g_parameterisation[
      self.scaler.active_parameterisation].calculate_smooth_derivatives()
    return derivatives'''
