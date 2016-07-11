#!/usr/bin/env cctbx.python
#
#  Copyright (C) (2015) Diamond Light Source
#
#  Author: Graeme Winter
#
#  This code is distributed under the BSD license, a copy of which is
#  included in the root directory of this package.
#

# python imports
from __future__ import division
import os
import libtbx.load_env # required for libtbx.env.find_in_repositories
from libtbx import easy_run
from libtbx.test_utils import open_tmp_directory

def test1():

  dials_regression = libtbx.env.find_in_repositories(
    relative_path="dials_regression",
    test=os.path.isdir)

  data_dir = os.path.join(dials_regression, "misc_test_data")
  experiments_path = os.path.join(data_dir, "i04-indexed.json")
  pickle_path = os.path.join(data_dir, "i04-indexed.pickle")

  for path in (experiments_path, pickle_path):
    assert os.path.exists(path)


  cmd = "dials.check_indexing_symmetry " + \
    "%s %s d_min=4 d_max=10 symop_threshold=0.6" % \
    (experiments_path, pickle_path)

  # work in a temporary directory

  cwd = os.path.abspath(os.curdir)
  tmp_dir = open_tmp_directory(suffix="test_dials_refine")
  os.chdir(tmp_dir)
  try:
    result = easy_run.fully_buffered(command=cmd).raise_if_errors()
  finally:
    os.chdir(cwd)

  print "OK"
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
