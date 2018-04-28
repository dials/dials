"""Tests for ObservationManager and related objects."""

from __future__ import absolute_import, division, print_function
from libtbx.test_utils import approx_equal

from cctbx.array_family import flex
from dials_scaling_helpers_ext import GroupedObservations, \
  minimum_multiplicity_selection
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


def test_grouped_observations():
  # some dummy observations for one reflection
  I1, I2, I3 = 101, 100, 99
  w1, w2, w3 = 0.9, 1.0, 0.8
  g1, g2, g3 = 1.01, 1.0, 0.99

  # combine with observations from some other group to make a dataset
  hkl = flex.miller_index([(0,0,1)] * 3 + [(0,0,2)] * 2)
  intensity = flex.double([I1, I2, I3, 10, 10])
  weight = flex.double([w1, w2, w3, 1.0, 1.0])
  phi = flex.double(len(hkl), 0) # dummy
  scale = flex.double([g1, g2, g3, 1.0, 1.0])

  # group the observations by Miller index
  go = GroupedObservations(hkl,
                           intensity,
                           weight,
                           phi,
                           scale)

  # ensure there are two groups, the first of size 3, the second of size 2
  assert list(go.get_group_size()) == [3,2]

  # the first group has an average intensity given by the HRS formula expanded:
  avI = (w1*g1*I1 + w2*g2*I2 + w3*g3*I3) / (w1*g1*g1 + w2*g2*g2 + w3*g3*g3)
  assert approx_equal(avI, go.get_average_intensity()[0])
