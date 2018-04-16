"""
Restraints manager classes for scaling.
"""
from scitbx import sparse
from dials.array_family import flex


class MultiScalingRestraints(object):
  """Scaling Restraints class for multi-dataset scaling."""

  def __init__(self, multi_apm):
    self.apm = multi_apm

  def calculate_restraints(self):
    """Calculate restraints for the scaling model."""
    R = flex.double([])
    G = flex.double([])
    for single_apm in self.apm.apm_list:
      restr = ScalingRestraints(single_apm).calculate_restraints()
      R.extend(restr[0])
      G.extend(restr[1])
    return [R, G]

  def calculate_jacobian_restraints(self):
    """Calculate restraints for jacobian."""
    residual_restraints = self.calculate_restraints()[0]
    n_restraints = residual_restraints.size() - residual_restraints.count(0.0)
    if n_restraints:
      weights = flex.double([])
      restraints_vector = flex.double([])
      jacobian = sparse.matrix(n_restraints, self.apm.n_active_params)
      cumul_restr_pos = 0
      for i, single_apm in enumerate(self.apm.apm_list):
        restraints = ScalingRestraints(single_apm).calculate_jacobian_restraints()
        if restraints:
          jacobian.assign_block(restraints[1], cumul_restr_pos, self.apm.apm_data[i]['start_idx'])
          cumul_restr_pos += restraints[1].n_rows
          restraints_vector.extend(restraints[0])
          weights.extend(restraints[2])
      return [restraints_vector, jacobian, weights]
    return None

class ScalingRestraints(object):
  """Scaling Restraints class."""

  def __init__(self, apm):
    self.apm = apm

  def calculate_jacobian_restraints(self):
    """Calculate restraints for jacobian."""
    residual_restraints = self.calculate_restraints()[0]
    n_restraints = residual_restraints.size() - residual_restraints.count(0.0)
    if n_restraints:
      weights = flex.double([])
      restraints_vector = flex.double([])
      jacobian = sparse.matrix(n_restraints, self.apm.n_active_params)
      cumul_restr_pos = 0
      for comp in self.apm.components.itervalues():
        restraints = comp['object'].calculate_jacobian_restraints()
        if restraints:
          jacobian.assign_block(restraints[1], cumul_restr_pos, comp['start_idx'])
          cumul_restr_pos += comp['n_params']
          restraints_vector.extend(restraints[0])
          weights.extend(restraints[2])
      # Return the restraints vector, jacobian and weights (unity as weights
      # contained in individual jacobian/weoghts calculations).
      return [restraints_vector, jacobian, weights]
    return None

  def calculate_restraints(self):
    """Calculate restraints for the scaling model."""
    residuals = flex.double([])
    gradient_vector = flex.double([])
    for comp in self.apm.components.itervalues():
      resid = comp['object'].calculate_restraints()
      if resid:
        gradient_vector.extend(resid[1])
        residuals.extend(resid[0])
      else:
        gradient_vector.extend(flex.double(comp['n_params'], 0.0))
    return [residuals, gradient_vector]
