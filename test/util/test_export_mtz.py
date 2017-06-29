"""
Unit testing for the export_mtz.py routines
"""

from __future__ import absolute_import, division

try:
  from mock import Mock
except ImportError:
  from unittest.mock import Mock

import pytest
import itertools

import dials.util.export_mtz as export_mtz

def in_ranges(value, ranges):
  """Check if a value is in a list of ranges"""
  return all( low <= value <= high for low, high in ranges)

def has_consecutive_ranges(ranges):
  """Check that a set of ranges is non-consecutive"""
  total_entries = sum(h-l+1+1 for l, h in ranges)
  union = set().union(*[set(range(l,h+2)) for l,h in ranges])
  return len(union) < total_entries

def offset_ranges(offsets, ranges):
  return [(l+off, h+off) for off, (l, h) in zip(offsets, ranges)]

class TestBatchRangeCalculations(object):
  "Test the calculation of non-overlapping batch ranges"

  class MockExperiment(object):
    def __init__(self, image_range):
      assert len(image_range) == 2
      self.image_range = tuple(image_range)

  def _run_ranges(self, ranges):
    """Convenience method to run the routine with a minimal experiment, and return the result as ranges"""
    input = [self.MockExperiment(x) for x in ranges]
    return offset_ranges(export_mtz._calculate_batch_offsets(input), ranges)

  def test_calculate_batch_ranges(self):
    assert self._run_ranges([(1,1)]) == [(1,1)]
    # Zero is shifted
    assert not in_ranges(0, self._run_ranges([(0,0)])), "Should be no zeroth batch"
    
    input = [self.MockExperiment(x) for x in [(1,1), (1,1)]]
    assert not set(self._run_ranges([(1,1), (1,1)])) == {(1,1)}, "Overlapping simple ranges"
    
    data_tests = [
      [(1,1), (1,1)],
      [(1,1), (8,8), (9,9)],
      [(23,24), (70, 100), (1,1), (1,4), (1,1)]
    ]
    for data in data_tests:
      print("Running ", data)
      print("  ", self._run_ranges(data))
      assert all([float(x).is_integer() for x in itertools.chain(*self._run_ranges(data))]), "Fractional epochs"
      assert all(isinstance(x, int) for x in itertools.chain(*self._run_ranges(data))), "Not all true integers"
      assert not has_consecutive_ranges(self._run_ranges(data))


