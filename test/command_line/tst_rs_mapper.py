#!/usr/bin/env cctbx.python

#
#  Copyright (C) (2014) STFC Rutherford Appleton Laboratory, UK.
#
#  Author: David Waterman.
#
#  This code is distributed under the BSD license, a copy of which is
#  included in the root directory of this package.
#

"""
Test dials.rs_mapper on images if dials_regression is present.

"""

# python imports
from __future__ import absolute_import, division

import os

import libtbx.load_env  # required for libtbx.env.find_in_repositories
from libtbx import easy_run
from libtbx.test_utils import approx_equal, open_tmp_directory

def test1():

  dials_regression = libtbx.env.find_in_repositories(
    relative_path="dials_regression",
    test=os.path.isdir)

  data_dir = os.path.join(dials_regression, "centroid_test_data")
  datablock_path = os.path.join(data_dir, "datablock.json")

  # work in a temporary directory
  cwd = os.path.abspath(os.curdir)
  tmp_dir = open_tmp_directory(suffix="tst_rs_mapper")
  os.chdir(tmp_dir)
  cmd = 'dials.rs_mapper ' + datablock_path + ' map_file="junk.ccp4"'

  result = easy_run.fully_buffered(command=cmd).raise_if_errors()
  # load results
  from iotbx import ccp4_map
  from scitbx.array_family import flex
  m = ccp4_map.map_reader(file_name="junk.ccp4")
  assert len(m.data) == 7189057
  assert approx_equal(m.header_min, -1.0)
  assert approx_equal(flex.min(m.data), -1.0)

  assert approx_equal(m.header_max, 2052.75)
  assert approx_equal(flex.max(m.data), 2052.75)

  assert approx_equal(m.header_mean, 0.018606403842568398)
  assert approx_equal(flex.mean(m.data), 0.018606403842568398)


  return

def run():
  if not libtbx.env.has_module("dials_regression"):
    print "Skipping tests in " + __file__ + " as dials_regression not present"
    return

  test1()

if __name__ == '__main__':
  from dials.test import cd_auto
  with cd_auto(__file__):
    from libtbx.utils import show_times_at_exit
    show_times_at_exit()
    run()
