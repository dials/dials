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

class ScaleFactor(object):
  """Base class for parameterisation of a component of the overall scale"""

  def get_param_vals(self):
    """export the values of the internal list of parameters as a
    sequence of floats."""

    # In this case, trivial as there is only one parameter set and it is not
    # fixed, but always free
    return self._param.value

  def set_param_vals(self, vals):
    """set the values of the internal list of parameters from a
    sequence of floats."""

    # In this case, trivial as there is only one parameter set and it is not
    # fixed, but always free
    self._param.value = vals

    return

  def get_factors_and_derivatives(self, seq):
    """Calculate and return the smoothed values and their derivatives with
    respect to the underlying parameters for this scale factor component,
    for the reflections described by the values in seq. For example, these
    may be the phi rotation angle values for a scale factor component that is
    a function of phi."""

    # The Python version of the smoother is used like this for a single smoothed
    # value:
    #
    # val, weights, sumweight = self._smoother.value_weight(phi, param_set)
    #
    # Where val and sumweight are both single values, val being the smoothed
    # value at position 'phi', and sumweight is the sum of weights. Weights
    # is a sparse vector typically with 3 non-zero values.
    #
    # The gradient of the smoothed value with respect to the parameters,
    # d[val]/d[p], is given by
    #
    # dvp_dp = weights * (1. / sumweight)
    #
    # We can start by looping over the reflection data (values in 'seq') and
    # performing these calculations, but it will be very slow. Ultimately it
    # would be better for the smoother to move to C++ and for these calculations
    # to be available in a 'vectorised' way so that from Python a single call
    # is made to the smoother and the loop over seq is in C++. In that case,
    # val and sumweight become vectors of length seq, while dv_dp is a matrix.
    # This matrix is sparse, so better to use the sparse matrix type. A decision
    # is required as to whether the weights for a single reflection span the
    # rows or columns of this matrix.
    #
    # I think that eventually the derivatives of the overall scale factor are
    # needed in such a way that each row is a single reflection, so arrange it
    # this way.

    pass

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

