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
from libtbx.test_utils import approx_equal

from scitbx.array_family import flex
from scitbx import sparse
from dials_scaling_helpers_ext import row_multiply
from dials.algorithms.scaling.scale_parameterisation import ScaleParameterisation

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

if __name__ == '__main__':

  test_row_multiply()
