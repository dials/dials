#!/usr/bin/env cctbx.python

#
#  Copyright (C) (2016) STFC Rutherford Appleton Laboratory, UK.
#
#  Author: David Waterman.
#
#  This code is distributed under the BSD license, a copy of which is
#  included in the root directory of this package.
#

"""
Test command line program dials.two_theta_refine by running a job with saved
data and comparing with expected output.

"""

# python imports
from __future__ import division
import os
import shutil
import libtbx.load_env # required for libtbx.env.find_in_repositories
from libtbx import easy_run
from libtbx.test_utils import open_tmp_directory, approx_equal
from dxtbx.model.experiment.experiment_list import ExperimentListFactory
from dials.array_family import flex

def test1():

  dials_regression = libtbx.env.find_in_repositories(
    relative_path="dials_regression",
    test=os.path.isdir)

  # use multiple scan small molecule data for this test
  data_dir = os.path.join(dials_regression, "xia2-28")
  prefix = ["20", "25", "30", "35"]
  exp_path = [e + "_integrated_experiments.json" for e in prefix]
  exp_path = [os.path.join(data_dir, e) for e in exp_path]
  pkl_path = [e + "_integrated.pickle" for e in prefix]
  pkl_path = [os.path.join(data_dir, e) for e in pkl_path]

  for pth in exp_path + pkl_path:
    assert os.path.exists(pth)

  cmd = "dials.two_theta_refine " + " ".join(exp_path) + " " + " ".join(pkl_path)
  print cmd

  # work in a temporary directory
  cwd = os.path.abspath(os.curdir)
  tmp_dir = open_tmp_directory(suffix="test_dials_two_theta_refine")
  os.chdir(tmp_dir)
  try:
    result = easy_run.fully_buffered(command=cmd).raise_if_errors()
    ref_exp = ExperimentListFactory.from_json_file("refined_cell.json",
                check_format=False)
  finally:
    os.chdir(cwd)
    # clean up tmp dir
    shutil.rmtree(tmp_dir)

  xls = ref_exp.crystals()
  assert len(xls) == 1 # crystal models should have been combined
  xl = xls[0]

  # test refined crystal model against expected values
  assert approx_equal(xl.get_unit_cell().parameters(),
    (5.428022880, 8.144145476, 12.039666971, 90.0, 90.0, 90.0))
  assert approx_equal(xl.get_cell_parameter_sd(),
    (9.58081e-5, 0.000149909, 0.000215765, 0, 0, 0))
  assert approx_equal(xl.get_cell_volume_sd(), 0.0116254298)

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
