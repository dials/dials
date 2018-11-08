"""
Test for the functions in dials.util.batch_handling.

These functions are designed to handle reflection/experiment lists that
include both sweep experiments and single image datasets, which do not
have a scan object.
"""

from dials.util.batch_handling import exclude_batches_in_reflections, \
  assign_batches_to_reflections, get_current_batch_ranges_for_scaling, \
  calculate_batch_offsets, set_batch_offsets, \
  get_batch_ranges, get_image_ranges, _next_epoch
from mock import Mock
from dials.array_family import flex
from dxtbx.model import Experiment, Scan
import pytest

def reflections_1():
  """Test reflection table with batch"""
  r = flex.reflection_table()
  r['batch'] = flex.int([1, 11, 21, 31, 41, 51, 61, 71, 81, 91])
  r.set_flags(flex.bool(10, False), r.flags.user_excluded_in_scaling)
  return r

def reflections_2():
  """Test reflection table with batch"""
  r = flex.reflection_table()
  r['batch'] = flex.int([201, 211, 221, 231, 241, 251, 261, 271, 281, 291])
  r.set_flags(flex.bool(10, False), r.flags.user_excluded_in_scaling)
  return r

def reflections_3():
  """Test reflection table with xyzobs.px.value"""
  r = flex.reflection_table()
  r['xyzobs.px.value'] = flex.vec3_double([(0, 0, 0.5), (0, 0, 1.5)])
  return r

def mock_experiments():
  """Mock experiments list with 2 experiments"""
  explist = []
  explist.append(Mock())
  explist[0].identifier = "0"
  explist.append(Mock())
  explist[1].identifier = "1"
  return explist

def mock_experiment_with_scaling_model(batch_range=(1, 50)):
  """Mock experiment with a scaling model valid_batch_range"""
  exp = Mock()
  exp.scaling_model.configdict = {'valid_batch_range': batch_range}
  return exp

def mock_experiment_without_scaling_model():
  """Mock experiment with no scaling model or scan"""
  exp = Mock()
  exp.scaling_model = None
  return exp

def test_assign_batches_to_reflections():
  """Test for namesake function"""
  reflections = [reflections_3(), reflections_3()]
  reflections = assign_batches_to_reflections(reflections, batch_offsets=[0, 100])
  assert list(reflections[0]['batch']) == [1, 2]
  assert list(reflections[1]['batch']) == [101, 102]

def test_exclude_batches_in_reflections():
  """Tests for namesake function"""
  # Simple case
  explist = mock_experiments()
  reflections = [reflections_1(), reflections_2()]
  in_use_batch_ranges = [(1, 100), (201, 300)]
  exclude_batches = [[251, 300]]
  reflections, new_ranges = exclude_batches_in_reflections(
    reflections, explist, in_use_batch_ranges, exclude_batches)
  assert new_ranges[0] == (1, 100)
  assert new_ranges[1] == (201, 250)
  assert list(reflections[0].get_flags(reflections[0].flags.user_excluded_in_scaling)) == \
    [False] * 10
  assert list(reflections[1].get_flags(reflections[1].flags.user_excluded_in_scaling)) == \
    [False] * 5 + [True] * 5

  # Case where excluding multiple ranges
  reflections = [reflections_1(), reflections_2()]
  exclude_batches = [[251, 300], [201, 220], [1, 20]]
  reflections, new_ranges = exclude_batches_in_reflections(
    reflections, explist, in_use_batch_ranges, exclude_batches)
  assert new_ranges[0] == (21, 100)
  assert new_ranges[1] == (221, 250)
  assert list(reflections[0].get_flags(reflections[0].flags.user_excluded_in_scaling)) == \
    [True] * 2 + [False] * 8
  assert list(reflections[1].get_flags(reflections[1].flags.user_excluded_in_scaling)) == \
    [True] * 2 + [False] * 3 + [True] * 5

  # Case where excluding in middle - returns initial range but excludes in middle
  reflections = [reflections_1(), reflections_2()]
  exclude_batches = [[251, 270]]
  reflections, new_ranges = exclude_batches_in_reflections(
    reflections, explist, in_use_batch_ranges, exclude_batches)
  assert new_ranges[0] == (1, 100)
  assert new_ranges[1] == (201, 300)
  assert list(reflections[0].get_flags(reflections[0].flags.user_excluded_in_scaling)) == \
    [False] * 10
  assert list(reflections[1].get_flags(reflections[1].flags.user_excluded_in_scaling)) == \
    [False] * 5 + [True] * 2 + [False] * 3

  with pytest.raises(ValueError):
    reflections, new_ranges = exclude_batches_in_reflections(
      reflections, explist, in_use_batch_ranges, [[1, 100]])

  #Test case where request outside of range - should not apply this exclusion
  exclude_batches = [[251, 350]]
  reflections = [reflections_1(), reflections_2()]
  reflections, new_ranges = exclude_batches_in_reflections(
    reflections, explist, in_use_batch_ranges, exclude_batches)
  assert new_ranges[0] == (1, 100)
  assert new_ranges[1] == (201, 300)
  assert list(reflections[0].get_flags(reflections[0].flags.user_excluded_in_scaling)) == \
    [False] * 10
  assert list(reflections[1].get_flags(reflections[1].flags.user_excluded_in_scaling)) == \
    [False] * 10

def test_get_current_batch_ranges_for_scaling():
  """Test for namesake function"""
  batch_ranges = [(1, 100), (1, 50), (1, 200)]
  exp = mock_experiment_with_scaling_model(batch_range=(1, 50))
  exp2 = mock_experiment_with_scaling_model(batch_range=(1, 25))
  exp3 = mock_experiment_without_scaling_model()
  in_use_ranges = get_current_batch_ranges_for_scaling([exp, exp2, exp3], batch_ranges)
  assert in_use_ranges == [(1, 50), (1, 25), (1, 200)]

def test_calculate_batch_offsets():
  """Test offset calculation. Offset is next number ending in 01 bigger than
  previous batch numbers which is not consecutive"""
  scan = Scan(image_range=[1, 200], oscillation=[0.0, 1.0])
  exp1 = Experiment(scan=scan)
  exp2 = Experiment()
  offsets = calculate_batch_offsets([exp1, exp2])
  assert offsets == [0, 301]

def test_set_batch_offsets():
  """Test for namesake function"""
  scan = Scan(image_range=[1, 200], oscillation=[0.0, 1.0])
  exp1 = Experiment(scan=scan)
  exp2 = Experiment()
  set_batch_offsets([exp2, exp1], [0, 101])
  assert exp1.scan.get_batch_offset() == 101

def test_get_batch_ranges():
  """Test for namesake function"""
  scan = Scan(image_range=[1, 200], oscillation=[0.0, 1.0])
  exp1 = Experiment(scan=scan)
  exp2 = Experiment(scan=scan)
  batch_offsets = [0, 300]
  experiments = [exp1, exp2]
  batch_ranges = get_batch_ranges(experiments, batch_offsets)
  assert batch_ranges == [(1, 200), (301, 500)]

def test_get_image_ranges():
  """Test for namesake function"""
  scan = Scan(image_range=[1, 200], oscillation=[0.0, 1.0])
  exp1 = Experiment(scan=scan)
  exp2 = Experiment()
  experiments = [exp1, exp2]
  image_ranges = get_image_ranges(experiments)
  assert image_ranges == [(1, 200), (0, 0)]

def test_next_epoch():
  """Test for namesake function"""
  assert _next_epoch(100) == 201
  assert _next_epoch(99) == 101
  assert _next_epoch(105) == 201
