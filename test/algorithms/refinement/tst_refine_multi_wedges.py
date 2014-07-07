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
Test multiple wedge refinement.

"""

# python imports
from __future__ import division
import os
import libtbx.load_env # required for libtbx.env.find_in_repositories
from libtbx import easy_run
from libtbx.test_utils import approx_equal

def test1():

  dials_regression = libtbx.env.find_in_repositories(
    relative_path="dials_regression",
    test=os.path.isdir)

  data_dir = os.path.join(dials_regression, "refinement_test_data",
                          "multi_narrow_wedges")
  phil_path = os.path.join(data_dir, "all.phil")
  assert os.path.exists(phil_path)

  cwd = os.path.abspath(os.curdir)
  os.chdir(data_dir)
  command = "dials.python refine_multi_wedges.py all.phil"
  try:
    result = easy_run.fully_buffered(command=command).raise_if_errors()
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
  import sys
  from libtbx.utils import show_times_at_exit
  show_times_at_exit()
  run()
