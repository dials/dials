from __future__ import absolute_import, division, print_function

import os

from libtbx import easy_run
from dials.array_family import flex

def test_filter_reflections(tmpdir):
  tmpdir.chdir()

  # make a dummy reflection table for the test, setting some flags
  rt = flex.reflection_table.empty_standard(6)
  rt['iobs'] = flex.size_t_range(len(rt))
  mask1 = flex.bool([True] * 3 + [False] * 3)
  mask2 = flex.bool([True, False] * 3)
  rt.set_flags(mask1, rt.flags.integrated)
  rt.set_flags(mask2, rt.flags.reference_spot)
  rt_name = "test_refs.pickle"
  rt.as_pickle(rt_name)

  cmd = ("dev.dials.filter_reflections " + rt_name + " flag_expression="
         "'integrated & ~reference_spot'")

  result = easy_run.fully_buffered(command=cmd).raise_if_errors()
  # load results
  ref = flex.reflection_table.from_pickle("filtered.pickle")

  # The test selects only 1 reflection
  assert len(ref) == 1
  assert list(ref['iobs']) == [1]
