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
Test command line program dials.slice_sweep by running a job with saved data and
comparing with expected output.

"""

# python imports
from __future__ import absolute_import, division
import os
import libtbx.load_env # required for libtbx.env.find_in_repositories
from libtbx import easy_run
from libtbx.test_utils import open_tmp_directory
from dxtbx.model.experiment_list import ExperimentListFactory
import cPickle as pickle

def test1():

  dials_regression = libtbx.env.find_in_repositories(
    relative_path="dials_regression",
    test=os.path.isdir)

  # use the i04_weak_data for this test
  data_dir = os.path.join(dials_regression, "refinement_test_data", "i04_weak_data")
  experiments_path = os.path.join(data_dir, "experiments.json")
  pickle_path = os.path.join(data_dir, "indexed_strong.pickle")

  for pth in (experiments_path, pickle_path):
    assert os.path.exists(pth)

  cmd = "dials.slice_sweep " + experiments_path + " " + pickle_path + \
  ' "scan_range=1 20"'
  print cmd

  # work in a temporary directory
  cwd = os.path.abspath(os.curdir)
  tmp_dir = open_tmp_directory(suffix="test_dials_slice_sweep")
  os.chdir(tmp_dir)
  try:
    result = easy_run.fully_buffered(command=cmd).raise_if_errors()
    # load results
    sliced_exp = ExperimentListFactory.from_json_file("experiments_1_20.json",
                check_format=False)[0]
    with open("indexed_strong_1_20.pickle", "r") as f:
      sliced_refs = pickle.load(f)
  finally:
    os.chdir(cwd)

  # simple test of results
  assert sliced_exp.scan.get_image_range() == (1, 20)
  assert len(sliced_refs) == 3670

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
