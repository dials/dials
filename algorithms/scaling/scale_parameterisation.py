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
from scitbx import sparse
from dials.algorithms.refinement.parameterisation.scan_varying_model_parameters \
        import GaussianSmoother, ScanVaryingParameterSet
from dials_scaling_helpers_ext import row_multiply
from dials.algorithms.scaling.scaling_helpers import products_omitting_one_item

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

  def __len__(self):
    # There is only one parameter set, so the total number of parameters is
    # equal to the set length
    return self._set_len

class BFactor(ScaleFactor):
  """Smoothly varying B-factor describing falloff with resolution as a function
  of rotation angle (and hence addressing bulk radiation damage to some extent)
  """

  pass

class ScaleParameterisation(object):
  """Parameterisation of the overall scale, combining various separate factors,
  such as an incident beam factor, a B-factor correction and potentially other
  factors."""

  def __init__(self, factors_list):
    """Initialise with a list of component factors for the overall scale"""

    self._factors = factors_list

    # Check there are free parameters to refine
    self._length = sum(len(f) for f in self._factors)
    if self._length == 0:
      raise RuntimeError("There are no free parameters for refinement")

  def __len__(self):
    return self._length

  def get_param_vals(self):
    """Return a concatenated list of parameters from each of the components
    in the global model"""

    global_p_list = []
    for f in self._factors:
      global_p_list.extend(f.get_param_vals())

    return global_p_list

  def set_param_vals(self, vals):
    """Set the parameter values of the contained models to the values in
    vals. This list must be of the same length as the result of get_param_vals
    and must contain the parameter values in the same order."""

    assert len(vals) == len(self)
    it = iter(vals)

    for f in self._factors:
      f.set_param_vals([it.next() for i in xrange(len(f))])

    return

  def scales_and_derivatives(self, phi):
    """Calculate the overall scale factor at each position in 'phi' and the
    derivatives of that scale factor wrt all parameters of the model"""

    # obtain data from all scale factor components
    data = [f.get_factors_and_derivatives(phi) for f in self._factors]

    # FIXME only using phi at the moment. In future will need other information
    # such as s1 directions or central impacts, and will have to pass the right
    # bits of information to the right ScaleFactor components contained here
    scale_components, grad_components = zip(*data)

    # the overall scale is the product of the separate scale components
    overall_scale = reduce(lambda fac1, fac2: fac1 * fac2, scale_components)

    # to convert derivatives of each scale component to derivatives of the
    # overall scale we multiply by the product of the other scale components,
    # omitting the factor that has been differentiated
    if len(scale_components) > 1:
      omit_one_prods = products_omitting_one_item(scale_components)

      grad_components = [row_multiply(g, coeff) for g, coeff in \
                         zip(grad_components, omit_one_prods)]

    # Now combine the gradient components by columns to produce a single
    # gradient matrix for the overall scale
    each_ncol = [g.n_cols for g in grad_components]
    tot_ncol = sum(each_ncol)
    grad = sparse.matrix(len(overall_scale), tot_ncol)
    col_start = [0] + each_ncol[:-1]
    for icol, g in zip(col_start, grad_components):
      grad.assign_block(g, 0, icol)

    return overall_scale, grad

