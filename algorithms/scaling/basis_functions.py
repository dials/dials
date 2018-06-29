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
  def __init__(self, curvatures=False):
    self.curvatures = curvatures

  def calc_component_scales_derivatives(self, apm, block_id):
    """Calculate the scales and derivatives for all components for a given
    block, returning each as a list of values from the components."""
    scales = []
    derivatives = []
    curvatures = []
    for component in apm.components.itervalues():
      sdc = component['object'].calculate_scales_and_derivatives(block_id, self.curvatures)
      scales.append(sdc[0])
      derivatives.append(sdc[1])
      if self.curvatures:
        curvatures.append(sdc[2])
    return scales, derivatives, curvatures

  @staticmethod
  def calculate_scale_factors(apm, block_id, scales):
    """Calculate the overall scale factor for each reflection from individual
    components."""
    #msf_list = []
    #for block_id in range(len(apm.n_obs)):
    if not scales:
      return None
    multiplied_scale_factors = flex.double(scales[0].size(), 1.0)
    for s in scales:
      multiplied_scale_factors *= s
    #for component in apm.components.itervalues():
    #  multiplied_scale_factors *= component['object'].inverse_scales[block_id]
    if apm.constant_g_values:
      multiplied_scale_factors *= apm.constant_g_values[block_id]
    return multiplied_scale_factors
    #msf_list.append(multiplied_scale_factors)
    #return msf_list

  @staticmethod
  def calculate_derivatives(apm, block_id, scales, derivatives_list):
    """Calculate the derivatives matrix."""
    if not scales:
      return None
    if len(scales) == 1:
      #only one active parameter, so don't need to chain rule any derivatives
      #for block_id in range(len(apm.n_obs)):
      return derivatives_list[0]
        #deriv_list.append(apm.components.values()[0]['object'].derivatives[block_id])
      #return deriv_list
    #for block_id in range(len(apm.n_obs)):
    derivatives = sparse.matrix(apm.n_obs[block_id], apm.n_active_params)
    col_idx = 0
    for i, d in enumerate(derivatives_list):
      scale_multipliers = flex.double(apm.n_obs[block_id], 1.0)
      for j, s1 in enumerate(scales):
        if i != j:
          scale_multipliers *= s1
      if apm.constant_g_values:
        scale_multipliers *= apm.constant_g_values[block_id]
      next_deriv = row_multiply(d, scale_multipliers)
      derivatives.assign_block(next_deriv, 0, col_idx)
      col_idx += d.n_cols
    return derivatives
    '''for comp_name, component in apm.components.iteritems():
      derivs = component['object'].derivatives[block_id]
      scale_multipliers = flex.double(apm.n_obs[block_id], 1.0)
      for comp, obj in apm.components.iteritems():
        if comp != comp_name:
          scale_multipliers *= obj['object'].inverse_scales[block_id]
      if apm.constant_g_values:
        scale_multipliers *= apm.constant_g_values[block_id]
      next_deriv = row_multiply(derivs, scale_multipliers)
      derivatives.assign_block(next_deriv, 0, component['start_idx'])
    return derivatives'''
    #deriv_list.append(derivatives)
    #return deriv_list

  @staticmethod
  def calculate_curvatures(apm, block_id, scales, curvatures_list):
    """Calculate the curvatures matrix."""
    if not scales:
      return None
    if len(scales) == 1:
      return curvatures_list[0]
    curvatures = sparse.matrix(apm.n_obs[block_id], apm.n_active_params)
    col_idx = 0
    for i, c in enumerate(curvatures_list):
      scale_multipliers = flex.double(apm.n_obs[block_id], 1.0)
      for j, s1 in enumerate(scales):
        if i != j:
          scale_multipliers *= s1
      if apm.constant_g_values:
        scale_multipliers *= apm.constant_g_values[block_id]
      next_curv = row_multiply(c, scale_multipliers)
      curvatures.assign_block(next_curv, 0, col_idx)
      col_idx += c.n_cols
    return curvatures

  def calculate_scales_and_derivatives(self, apm, block_id):
    """Calculate and return scale factors, derivatives and optionally
    curvatures to be used in minimisation."""
    scales, derivatives, curvs = self.calc_component_scales_derivatives(apm, block_id)
    #self.update_scale_factors(apm, block_id)
    if self.curvatures:
      return (self.calculate_scale_factors(apm, block_id, scales),
       self.calculate_derivatives(apm, block_id, scales, derivatives),
       self.calculate_curvatures(apm, block_id, scales, curvs))
    return self.calculate_scale_factors(apm, block_id, scales), \
      self.calculate_derivatives(apm, block_id, scales, derivatives)
