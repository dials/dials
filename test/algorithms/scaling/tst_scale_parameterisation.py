#!/usr/bin/env cctbx.python

#
#  Copyright (C) (2016) STFC Rutherford Appleton Laboratory, UK.
#
#  Author: David Waterman.
#
#  This code is distributed under the BSD license, a copy of which is
#  included in the root directory of this package.
#

"""Tests for ScaleParameterisation and related objects."""

from __future__ import division
import random

from scitbx.array_family import flex
from scitbx import sparse
from libtbx.test_utils import approx_equal

from dials_scaling_helpers_ext import row_multiply
from dials.algorithms.scaling.scale_parameterisation import ScaleParameterisation
from dials.algorithms.scaling.scale_parameterisation import IncidentBeamFactor

def test_row_multiply():

  m = sparse.matrix(3, 2)
  m[0,0] = 1.
  m[0,1] = 2.
  m[1,1] = 3.
  m[2,0] = 4.

  fac = flex.double((3, 2, 1))

  m2 = row_multiply(m, fac)

  assert m2.as_dense_matrix().as_1d() == flex.double(
    [3.0, 6.0, 0.0, 6.0, 4.0, 0.0]).all_eq(True)
  print "OK"

class TestIncidentBeamFactor(object):

  def __init__(self):

    self.ibf = IncidentBeamFactor([0,180])
    assert self.ibf.get_param_vals() == [1] * 38

    # set the smoother parameters to something other than unity throughout
    v2 = [random.uniform(0.8, 1.2) for v in self.ibf.get_param_vals()]
    self.ibf.set_param_vals(v2)

    # request values at 3 phi positions, one of them randomly chosen
    self.phi = [0, random.uniform(0, 180), 180]

    v, dv_dp = self.ibf.get_factors_and_derivatives(self.phi)

    # calculate finite diff gradients. This comes back as a list of flex.doubles
    fd_dv_dp = self._calc_fd_grad()

    # convert to list of lists and construct 2D matrix
    fd_dv_dp = flex.double([list(e) for e in fd_dv_dp]).matrix_transpose()

    # compare with the analytical calculation, converted to dense
    assert approx_equal(dv_dp.as_dense_matrix(), fd_dv_dp)

    print "OK"

  def _calc_fd_grad(self):
    delta = 1.e-7
    fd_grad = []
    p_vals = self.ibf.get_param_vals()

    for i, p in enumerate(p_vals):
      p_vals[i] = p - delta/2
      self.ibf.set_param_vals(p_vals)
      rev_state = self.ibf.get_factors_and_derivatives(self.phi)[0]

      p_vals[i] = p + delta/2
      self.ibf.set_param_vals(p_vals)
      fwd_state = self.ibf.get_factors_and_derivatives(self.phi)[0]

      p_vals[i] = p
      self.ibf.set_param_vals(p_vals)

      fd_grad.append((fwd_state - rev_state)/delta)

    return fd_grad

def test_scale_parameterisation():

  ibf = IncidentBeamFactor([0,180])
  sf = ScaleParameterisation(factors_list=[ibf])

  # test getting and setting parameter values
  p = sf.get_param_vals()
  p2 = [random.uniform(0.8, 1.2) * e for e in p]
  sf.set_param_vals(p2)
  p3 = sf.get_param_vals()
  for e1, e2 in zip(p2, p3): assert e1 == e2
  print "OK"

  # test getting overall scale and its derivatives
  phi = [0, random.uniform(0, 180), 180]
  a, b = sf.scales_and_derivatives(phi)

  print "OK"
  return


if __name__ == '__main__':

  test_row_multiply()

  TestIncidentBeamFactor()

  test_scale_parameterisation()
