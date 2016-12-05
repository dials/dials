#!/usr/bin/env cctbx.python
#
#  Copyright (C) (2016) STFC Rutherford Appleton Laboratory, UK.
#
#  Author: David Waterman.
#
#  This code is distributed under the BSD license, a copy of which is
#  included in the root directory of this package.
#

"""Contains classes required for parameterising the scale factor components
and the overall scale."""

from __future__ import division
from dials.array_family import flex
from libtbx import phil
from libtbx.utils import Sorry
from dials.algorithms.refinement.parameterisation.scan_varying_model_parameters \
        import GaussianSmoother, ScanVaryingParameterSet
from dials_scaling_helpers_ext import row_multiply

class ScaleFactor(object):
  """Base class for parameterisation of a component of the overall scale"""

  def get_param_vals(self):
    """export the values of the internal list of parameters as a
    sequence of floats."""

    # In this case, trivial as there is only one parameter set and it is not
    # fixed, but always free. Make sure there is a copy
    return list(self._param.value)

  def set_param_vals(self, vals):
    """set the values of the internal list of parameters from a
    sequence of floats."""

    # In this case, trivial as there is only one parameter set and it is not
    # fixed, but always free. Make sure there is a copy
    self._param.value = list(vals)

    return

  def get_factors_and_derivatives(self, seq):
    """Calculate and return the smoothed values and their derivatives with
    respect to the underlying parameters for this scale factor component,
    for the reflections described by the values in seq. For example, these
    may be the phi rotation angle values for a scale factor component that is
    a function of phi."""

    # Obtain data from the smoother, where value and sumweight are arrays with
    # the same length as seq. Weight is a sparse matrix with one row per
    # element of seq. Each row has some (typically 3) non-zero values
    # corresponding to the positions nearest the element of seq.
    value, weight, sumweight = self._smoother.multi_value_weight(seq,
      self._param)

    # The gradient of a single smoothed value with respect to the parameters,
    # d[v]/d[p] = w_i * (1. / sum(w_i)), where the sum is over the parameters.
    # Calculate this for all values v.
    inv_sw = 1. / sumweight
    dv_dp = row_multiply(weight, inv_sw)

    return value, dv_dp

class IncidentBeamFactor(ScaleFactor):
  """Smoothly varying scale factor combining incident beam flux variation
  with other factors that can be modelled as a function of rotation angle"""

  def __init__(self, phi_range_deg, deg_per_interval=5):

    n_intervals = max(int(abs(
      phi_range_deg[1] - phi_range_deg[0]) / deg_per_interval), 1)

    # Set up the smoother
    self._smoother = GaussianSmoother(phi_range_deg, n_intervals)
    self._set_len = self._smoother.num_values()

    # initial value of scale factors is 1
    value = 1
    self._param = ScanVaryingParameterSet(value, self._set_len,
      name = "IncidentBeam")

class BFactor(ScaleFactor):
  """Smoothly varying B-factor describing falloff with resolution as a function
  of rotation angle (and hence addressing bulk radiation damage to some extent)
  """

  pass

class ScaleParameterisation(object):
  """Parameterisation of the overall scale, combining an incident beam factor
  with a B-factor correction and potentially other factors."""
  pass

if __name__ == '__main__':

  ibf = IncidentBeamFactor([0,180])
  assert ibf.get_param_vals() == [1] * 38

