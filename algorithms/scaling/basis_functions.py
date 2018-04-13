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
  """Class that takes in a scaling_apm and calcuates the scale factors,
  derivatives and optionally curvatures for minimisation."""
  def __init__(self, apm, curvatures=False):
    self.apm = apm
    self.curvatures = curvatures

  def update_scale_factors(self):
    """Update the parameters in each SF object from the apm parameter list."""
    for component in self.apm.components.itervalues():
      component['object'].calculate_scales_and_derivatives(curvatures=self.curvatures)

  def calculate_scale_factors(self):
    """Calculate the overall scale factor for each reflection from individual
    components."""
    multiplied_scale_factors = flex.double(self.apm.n_obs, 1.0)
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
    derivatives = sparse.matrix(self.apm.n_obs, self.apm.n_active_params)
    for comp_name, component in self.apm.components.iteritems():
      derivs = component['object'].derivatives
      scale_multipliers = flex.double(self.apm.n_obs, 1.0)
      for comp, obj in self.apm.components.iteritems():
        if comp != comp_name:
          scale_multipliers *= obj['object'].inverse_scales
      if self.apm.constant_g_values:
        scale_multipliers *= self.apm.constant_g_values
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
    curvatures = sparse.matrix(self.apm.n_obs, self.apm.n_active_params)
    for comp_name, component in self.apm.components.iteritems():
      curvs = component['object'].curvatures
      if curvs != 0.0:
        scale_multipliers = flex.double(self.apm.n_obs, 1.0)
        for comp, obj in self.apm.components.iteritems():
          if comp != comp_name:
            scale_multipliers *= obj['object'].inverse_scales
        if self.apm.constant_g_values:
          scale_multipliers *= self.apm.constant_g_values
        next_curv = row_multiply(curvs, scale_multipliers)
        curvatures.assign_block(next_curv, 0, component['start_idx'])
    return curvatures

  def return_basis(self):
    """Calculate and return scale factors, derivatives and optionally
    curvatures to be used in minimisation."""
    self.update_scale_factors()
    if self.curvatures:
      return (self.calculate_scale_factors(), self.calculate_derivatives(),
        self.calculate_curvatures())
    return self.calculate_scale_factors(), self.calculate_derivatives()
