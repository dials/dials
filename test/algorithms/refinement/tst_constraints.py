#!/usr/bin/env cctbx.python

#
#  Copyright (C) (2017) STFC Rutherford Appleton Laboratory, UK.
#
#  Author: David Waterman.
#
#  This code is distributed under the BSD license, a copy of which is
#  included in the root directory of this package.
#

"""
Tests for the constraints system used in refinement

"""

from __future__ import absolute_import, division
from dials.array_family import flex
from dials.algorithms.refinement.constraints import EqualShiftConstraint
from dials.algorithms.refinement.constraints import ConstraintManager

def test1():

  x = flex.random_double(10)

  # constrain parameters 2 and 4 and 6, 7 and 8
  c1 = EqualShiftConstraint([1, 3], x)
  c2 = EqualShiftConstraint([5, 6, 7], x)

  cm = ConstraintManager([c1, c2], len(x))
  constrained_x = cm.constrain_parameters(x)

  # check the constrained parameters are as expected
  assert len(constrained_x) == 7
  assert constrained_x[5] == flex.mean(x.select([1, 3]))
  assert constrained_x[6] == flex.mean(x[5:8])

  # minimiser would modify the constrained parameters
  mod_constrained_x = constrained_x + 10.

  # check the expanded parameters are as expected
  expanded = cm.expand_parameters(mod_constrained_x)
  assert x + 10. == expanded

if __name__ == '__main__':

  # simple test of constraint manager
  test1()
