from __future__ import absolute_import, division, print_function

import os
import procrunner

def test_sort_intensities(dials_regression, run_in_tmpdir):
  result = procrunner.run_process([
      'dev.dials.sort_reflections',
      os.path.join(dials_regression, "centroid_test_data", "integrated.pickle"),
      'key=intensity.sum.value',
      'output=sorted1.pickle',
  ])
  assert result['exitcode'] == 0
  assert result['stderr'] == ''
  assert os.path.exists("sorted1.pickle")

  from dials.array_family import flex
  data = flex.reflection_table.from_pickle("sorted1.pickle")
  assert_sorted(data['intensity.sum.value'])

def test_reverse_sort_intensities(dials_regression, run_in_tmpdir):
  result = procrunner.run_process([
      'dev.dials.sort_reflections',
      os.path.join(dials_regression, "centroid_test_data", "integrated.pickle"),
      'output=sorted2.pickle',
      'key=intensity.sum.value',
      'reverse=True'
  ])
  assert result['exitcode'] == 0
  assert result['stderr'] == ''
  assert os.path.exists("sorted2.pickle")

  from dials.array_family import flex
  data = flex.reflection_table.from_pickle("sorted2.pickle")
  assert_sorted(data['intensity.sum.value'], reverse=True)

def test_default_sort_on_miller_index(dials_regression, run_in_tmpdir):
  result = procrunner.run_process([
      'dev.dials.sort_reflections',
      os.path.join(dials_regression, "centroid_test_data", "integrated.pickle"),
      'output=sorted3.pickle'
  ])
  assert result['exitcode'] == 0
  assert result['stderr'] == ''
  assert os.path.exists("sorted3.pickle")

  from dials.array_family import flex
  data = flex.reflection_table.from_pickle("sorted3.pickle")
  mi1 = data['miller_index']
  orig = flex.reflection_table.from_pickle(
      os.path.join(dials_regression, "centroid_test_data", "integrated.pickle")
  )
  mi2 = flex.miller_index(sorted(orig['miller_index']))
  assert mi1.all_eq(mi2)

def assert_sorted(data, reverse=False):
  assert len(data) > 0
  elem = data[0]
  for x in data:
    if reverse is True:
      assert x <= elem
    else:
      assert x >= elem
    elem = x
