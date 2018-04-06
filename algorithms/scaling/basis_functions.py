"""
Classes that take in a scaler and minimisation parameters and
return the scale factors and derivatives of the scale factors w.r.t.
the parameters
"""
from __future__ import print_function
from dials.array_family import flex
from dials_scaling_helpers_ext import row_multiply
from scitbx import sparse

class basis_function(object):
  """Class that takes in a scaler and minimisation parameters and
  calcuates the scale factors, derivatives and optionally curvatures
  for minimisation."""
  def __init__(self, scaler, apm, curvatures=False):
    self.scaler = scaler
    self.apm = apm
    self.curvatures = curvatures

  def update_scale_factors(self):
    """Update the parameters in each SF object from the apm parameter list."""
    for component in self.apm.components:
      SF_object = self.apm.components[component]['object']
      SF_object.parameters = self.apm.select_parameters(component)
      SF_object.calculate_scales_and_derivatives()

  def update_scale_factors_with_curvs(self):
    """Update the parameters in each SF object from the apm parameter list."""
    for component in self.apm.components:
      SF_object = self.apm.components[component]['object']
      SF_object.parameters = self.apm.select_parameters(component)
      SF_object.calculate_scales_derivatives_curvatures()

  def calculate_scale_factors(self):
    """Calculate the overall scale factor for each reflection from individual
    components."""
    multiplied_scale_factors = flex.double(self.scaler.Ih_table.size, 1.0)
    for component in self.apm.components.itervalues():
      multiplied_scale_factors *= component['object'].inverse_scales
    if self.apm.constant_g_values:
      multiplied_scale_factors *= self.apm.constant_g_values
    return multiplied_scale_factors

  def calculate_derivatives(self):
    """Calculate the derivatives matrix."""
    if not self.apm.components:
      return None
    if len(self.apm.components) == 1:
      #only one active parameter, so don't need to chain rule any derivatives
      return self.apm.components.values()[0]['object'].derivatives
    derivatives = sparse.matrix(self.scaler.Ih_table.size,
      self.apm.n_active_params)
    for comp_name, component in self.apm.components.iteritems():
      derivs = component['object'].derivatives
      scale_multipliers = flex.double(self.scaler.Ih_table.size, 1.0)
      for comp, SF_obj in self.scaler.components.iteritems():
        if comp != comp_name:
          scale_multipliers *= SF_obj.inverse_scales
      next_deriv = row_multiply(derivs, scale_multipliers)
      derivatives.assign_block(next_deriv, 0, component['start_idx'])
    return derivatives

  def calculate_curvatures(self):
    """Calculate the curvatures matrix."""
    if not self.apm.components:
      return None
    if len(self.apm.components) == 1:
      #only one active parameter, so don't need to chain rule any derivatives
      return self.apm.components.values()[0]['object'].curvatures
    curvatures = sparse.matrix(self.scaler.Ih_table.size,
      self.apm.n_active_params)
    for comp_name, component in self.apm.components.iteritems():
      curvs = component['object'].curvatures
      if curvs != 0.0:
        scale_multipliers = flex.double(self.scaler.Ih_table.size, 1.0)
        for comp, SF_obj in self.scaler.components.iteritems():
          if comp != comp_name:
            scale_multipliers *= SF_obj.inverse_scales
        next_curv = row_multiply(curvs, scale_multipliers)
        curvatures.assign_block(next_curv, 0, component['start_idx'])
    return curvatures

  def return_basis(self):
    """Calculate and return scale factors, derivatives and optionally
    curvatures to be used in minimisation."""
    if self.curvatures:
      self.update_scale_factors_with_curvs()
      return (self.calculate_scale_factors(), self.calculate_derivatives(),
        self.calculate_curvatures())
    self.update_scale_factors()
    return self.calculate_scale_factors(), self.calculate_derivatives()

'''class xds_basis_function_log(basis_function):
  'Subclass of basis_function for logarithmic parameterisation'
  def calculate_scale_factors(self):
    SF_object = self.scaler.g_parameterisation[
      self.scaler.active_parameterisation]
    SF_object.set_scale_factors(self.apm.x)
    gxvalues = SF_object.calculate_smooth_scales()
    scale_factors = gxvalues + self.scaler.constant_g_values
    return flex.double(np.exp(scale_factors))

  def calculate_derivatives(self):
    'Derivatives are fixed by the parameterisation'
    derivatives = self.scaler.g_parameterisation[
      self.scaler.active_parameterisation].calculate_smooth_derivatives()
    return derivatives'''
