#!/usr/bin/env dials.python

#
#  Copyright (C) (2015) STFC Rutherford Appleton Laboratory, UK.
#
#  Author: David Waterman.
#
#  This code is distributed under the BSD license, a copy of which is
#  included in the root directory of this package.
#

"""
Test command line program dev.dials.filter_reflections by filtering according to
certain flags

"""

# python imports
from __future__ import division
import os
import shutil
from libtbx import easy_run
from libtbx.test_utils import open_tmp_directory
import cPickle as pickle

def test1():

  # work in a temporary directory
  cwd = os.path.abspath(os.curdir)
  tmp_dir = open_tmp_directory(suffix="test_dials_filter_reflections")
  os.chdir(tmp_dir)

  # make a dummy reflection table for the test, setting some flags
  from dials.array_family import flex
  rt = flex.reflection_table.empty_standard(6)
  rt['iobs'] = flex.size_t_range(len(rt))
  mask1 = flex.bool([True] * 3 + [False] * 3)
  mask2 = flex.bool([True, False] * 3)
  rt.set_flags(mask1, rt.flags.integrated)
  rt.set_flags(mask2, rt.flags.bad_spot)
  rt_name = "test_refs.pickle"
  rt.as_pickle(rt_name)

  cmd = "dev.dials.filter_reflections " + rt_name + " inclusions.flag=integrated" + \
    " exclusions.flag=bad_spot"
  print cmd

  try:
    result = easy_run.fully_buffered(command=cmd).raise_if_errors()
    # load results
    ref = flex.reflection_table.from_pickle("filtered.pickle")
  finally:
    os.chdir(cwd)

  # The test selects only 1 reflection
  assert len(ref) == 1
  assert list(ref['iobs']) == [1]

  print "OK"

  return

def run():

  test1()

if __name__ == '__main__':
  from dials.test import cd_auto
  with cd_auto(__file__):
    from libtbx.utils import show_times_at_exit
    show_times_at_exit()
    run()
