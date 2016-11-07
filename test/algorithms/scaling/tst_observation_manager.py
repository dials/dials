#!/usr/bin/env cctbx.python

#
#  Copyright (C) (2016) STFC Rutherford Appleton Laboratory, UK.
#
#  Author: David Waterman.
#
#  This code is distributed under the BSD license, a copy of which is
#  included in the root directory of this package.
#

"""Tests for ObservationManager and related objects."""

from __future__ import division
#from libtbx.test_utils import approx_equal

from cctbx.array_family import flex
from dials_scaling_helpers_ext import (GroupedObservations,
  minimum_multiplicity_selection)
from dials.algorithms.scaling.observation_manager import ObservationManager

def test_minimum_multiplicity_selection():

  # some groups of indices with increasing multiplicity
  hkl = flex.miller_index(
        [(0,0,1),
         (0,0,2), (0,0,2),
         (0,0,3), (0,0,3), (0,0,3),
         (0,0,4), (0,0,4), (0,0,4), (0,0,4),
         (0,0,5), (0,0,5), (0,0,5), (0,0,5), (0,0,5)])

  # test the various possible selections with multiplicity from 1 to 6
  sel = minimum_multiplicity_selection(hkl, 1)
  assert sel.count(True) == len(hkl)

  sel = minimum_multiplicity_selection(hkl, 2)
  assert (sel == flex.bool([False,] + [True] * 14)).all_eq(True)

  sel = minimum_multiplicity_selection(hkl, 3)
  assert (sel == flex.bool([False] * 3 + [True] * 12)).all_eq(True)

  sel = minimum_multiplicity_selection(hkl, 4)
  assert (sel == flex.bool([False] * 6 + [True] * 9)).all_eq(True)

  sel = minimum_multiplicity_selection(hkl, 5)
  assert (sel == flex.bool([False] * 10 + [True] * 5)).all_eq(True)

  sel = minimum_multiplicity_selection(hkl, 6)
  assert sel.count(False) == len(hkl)

  print "OK"

if __name__ == '__main__':

  test_minimum_multiplicity_selection()
