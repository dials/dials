#!/usr/bin/env cctbx.python

#
#  Copyright (C) (2013) STFC Rutherford Appleton Laboratory, UK.
#
#  Author: David Waterman.
#
#  This code is distributed under the BSD license, a copy of which is
#  included in the root directory of this package.
#

"""
Test command line program dials.refine by running a job with saved data and
comparing with expected output.

This serves as a high level test that not only checks whether refinement works,
but also that the command line program is functioning and that the output models
have not changed format and so on.

"""

# python imports
from __future__ import division
import os
import libtbx.load_env # required for libtbx.env.find_in_repositories
from libtbx import easy_run
from libtbx.test_utils import approx_equal
from libtbx.test_utils import open_tmp_directory

def test1():

  dials_regression = libtbx.env.find_in_repositories(
    relative_path="dials_regression",
    test=os.path.isdir)

  # use the i04_weak_data for this test
  data_dir = os.path.join(dials_regression, "refinement_test_data", "i04_weak_data")
  experiments_path = os.path.join(data_dir, "experiments.json")
  pickle_path = os.path.join(data_dir, "indexed_strong.pickle")

  for pth in ():
    assert os.path.exists(pth)

  cmd = "dials.refine " + experiments_path + " " + pickle_path

  # work in a temporary directory
  cwd = os.path.abspath(os.curdir)
  tmp_dir = open_tmp_directory(suffix="test_dials_refine")
  os.chdir(tmp_dir)
  result = easy_run.fully_buffered(command=cmd).raise_if_errors()
  os.chdir(cwd)

  # load results
  from dials.model.experiment import ExperimentListFactory
  reg_exp = ExperimentListFactory.from_json(
              os.path.join(data_dir, "regression_experiment.json"))
  ref_exp = ExperimentListFactory.from_json(
              os.path.join(data_dir, "refined_experiment.json"))

  # test refined models against expected
  assert reg_exp.crystal == ref_exp.crystal
  assert sweep.get_detector() == reg_sweep.get_detector()
  assert sweep.get_beam() == reg_sweep.get_beam()

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
