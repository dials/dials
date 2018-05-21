"""
Classes that take in a scaler and minimisation parameters and
return the scale factors and derivatives of the scale factors w.r.t.
the parameters
"""
from __future__ import print_function
from dials.array_family import flex
from dials_scaling_ext import row_multiply
from scitbx import sparse

class basis_function(object):
  """Class that takes in a scaling_apm and calcuates the scale factors,
  derivatives and optionally curvatures for minimisation."""
  def __init__(self, apm, curvatures=False):
    self.curvatures = curvatures
    self.apm = apm

  def update_scale_factors(self):
    """Update the parameters in each SF object from the apm parameter list."""
    for component in self.apm.components.itervalues():
      component['object'].calculate_scales_and_derivatives(self.curvatures)

  def calculate_scale_factors(self):
    """Calculate the overall scale factor for each reflection from individual
    components."""
    msf_list = []
    for block_id in range(len(self.apm.n_obs)):
      multiplied_scale_factors = flex.double(self.apm.n_obs[block_id], 1.0)
      for component in self.apm.components.itervalues():
        multiplied_scale_factors *= component['object'].inverse_scales[block_id]
      if self.apm.constant_g_values:
        multiplied_scale_factors *= self.apm.constant_g_values[block_id]
      msf_list.append(multiplied_scale_factors)
    return msf_list

  def calculate_derivatives(self):
    """Calculate the derivatives matrix."""
    if not self.apm.components:
      return None
    deriv_list = []
    if len(self.apm.components) == 1:
      #only one active parameter, so don't need to chain rule any derivatives
      for block_id in range(len(self.apm.n_obs)):
        deriv_list.append(self.apm.components.values()[0]['object'].derivatives[block_id])
      return deriv_list
    for block_id in range(len(self.apm.n_obs)):
      derivatives = sparse.matrix(self.apm.n_obs[block_id], self.apm.n_active_params)
      for comp_name, component in self.apm.components.iteritems():
        derivs = component['object'].derivatives[block_id]
        scale_multipliers = flex.double(self.apm.n_obs[block_id], 1.0)
        for comp, obj in self.apm.components.iteritems():
          if comp != comp_name:
            scale_multipliers *= obj['object'].inverse_scales[block_id]
        if self.apm.constant_g_values:
          scale_multipliers *= self.apm.constant_g_values[block_id]
        next_deriv = row_multiply(derivs, scale_multipliers)
        derivatives.assign_block(next_deriv, 0, component['start_idx'])
      deriv_list.append(derivatives)
    return deriv_list

  def calculate_curvatures(self):
    """Calculate the curvatures matrix."""
    if not self.apm.components:
      return None
    curv_list = []
    if len(self.apm.components) == 1:
      #only one active parameter, so don't need to chain rule any derivatives
      for block_id in range(len(self.apm.n_obs)):
        curv_list.append(self.apm.components.values()[0]['object'].curvatures[block_id])
      return curv_list
    for block_id in range(len(self.apm.n_obs)):
      curvatures = sparse.matrix(self.apm.n_obs[block_id], self.apm.n_active_params)
      for comp_name, component in self.apm.components.iteritems():
        curvs = component['object'].curvatures[block_id]
        if curvs != 0.0:
          scale_multipliers = flex.double(self.apm.n_obs[block_id], 1.0)
          for comp, obj in self.apm.components.iteritems():
            if comp != comp_name:
              scale_multipliers *= obj['object'].inverse_scales[block_id]
          if self.apm.constant_g_values:
            scale_multipliers *= self.apm.constant_g_values[block_id]
          next_curv = row_multiply(curvs, scale_multipliers)
          curvatures.assign_block(next_curv, 0, component['start_idx'])
      curv_list.append(curvatures)
    return curv_list

  def calculate_scales_and_derivatives(self):
    """Calculate and return scale factors, derivatives and optionally
    curvatures to be used in minimisation."""
    self.update_scale_factors()
    if self.curvatures:
      return (self.calculate_scale_factors(), self.calculate_derivatives(), None)
      #  self.calculate_curvatures(apm))
    return self.calculate_scale_factors(), self.calculate_derivatives()
