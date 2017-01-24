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
Test multiple stills refinement.

"""

# python imports
from __future__ import absolute_import, division
import os
import libtbx.load_env # required for libtbx.env.find_in_repositories
from libtbx import easy_run
from libtbx.test_utils import open_tmp_directory
from scitbx import matrix
from dxtbx.model.experiment.experiment_list import ExperimentListFactory

def test1():

  dials_regression = libtbx.env.find_in_repositories(
    relative_path="dials_regression",
    test=os.path.isdir)

  data_dir = os.path.join(dials_regression, "refinement_test_data",
                          "multi_stills")
  experiments_path = os.path.join(data_dir, "combined_experiments.json")
  reflections_path = os.path.join(data_dir, "combined_reflections.pickle")
  cmd = "dials.refine " + experiments_path + " " + reflections_path
  print cmd

  # work in a temporary directory
  cwd = os.path.abspath(os.curdir)
  tmp_dir = open_tmp_directory(suffix="tst_refine_multi_stills1")
  os.chdir(tmp_dir)
  try:
    result = easy_run.fully_buffered(command=cmd).raise_if_errors()
    # load results
    reg_exp = ExperimentListFactory.from_json_file(
                os.path.join(data_dir, "regression_experiments.json"),
                check_format=False)
    ref_exp = ExperimentListFactory.from_json_file("refined_experiments.json",
                check_format=False)
  finally:
    os.chdir(cwd)
  print "OK"

  # compare results
  tol = 1e-5
  for b1, b2 in zip(reg_exp.beams(), ref_exp.beams()):
    assert b1.is_similar_to(b2, wavelength_tolerance=tol,
                                direction_tolerance=tol,
                                polarization_normal_tolerance=tol,
                                polarization_fraction_tolerance=tol)
    s0_1 = matrix.col(b1.get_unit_s0())
    s0_2 = matrix.col(b2.get_unit_s0())
    assert s0_1.accute_angle(s0_2, deg=True) < 0.0057 # ~0.1 mrad
  for c1, c2 in zip(reg_exp.crystals(), ref_exp.crystals()):
    assert c1.is_similar_to(c2)

  for d1, d2 in zip(reg_exp.detectors(), ref_exp.detectors()):
    assert d1.is_similar_to(d2,
      fast_axis_tolerance=1e-4, slow_axis_tolerance=1e-4, origin_tolerance=1e-2)

  print "OK"

  return

def test2():
  """Compare results of multiprocess vs single process refinement to ensure
  they are the same"""

  dials_regression = libtbx.env.find_in_repositories(
    relative_path="dials_regression",
    test=os.path.isdir)

  data_dir = os.path.join(dials_regression, "refinement_test_data",
                          "multi_stills")
  experiments_path = os.path.join(data_dir, "combined_experiments.json")
  reflections_path = os.path.join(data_dir, "combined_reflections.pickle")
  cmd = "dials.refine " + experiments_path + " " + reflections_path + \
        " outlier.algorithm=null engine=LBFGScurvs output.reflections=None "
  cmd1 = cmd + "output.experiments=refined_experiments_nproc1.json nproc=1"
  print cmd1

  cmd2= cmd + "output.experiments=refined_experiments_nproc4.json nproc=4"
  print cmd2
  # work in a temporary directory
  cwd = os.path.abspath(os.curdir)
  tmp_dir = open_tmp_directory(suffix="tst_refine_multi_stills2")
  os.chdir(tmp_dir)
  try:
    result1 = easy_run.fully_buffered(command=cmd1).raise_if_errors()
    result2 = easy_run.fully_buffered(command=cmd2).raise_if_errors()
    # load results
    nproc1 = ExperimentListFactory.from_json_file(
      "refined_experiments_nproc1.json", check_format=False)
    nproc4 = ExperimentListFactory.from_json_file(
      "refined_experiments_nproc4.json", check_format=False)
  finally:
    os.chdir(cwd)
  print "OK"

  # compare results
  for b1, b2 in zip(nproc1.beams(), nproc4.beams()):
    assert b1.is_similar_to(b2)
  for c1, c2 in zip(nproc1.crystals(), nproc4.crystals()):
    assert c1.is_similar_to(c2)
  for d1, d2 in zip(nproc1.detectors(), nproc4.detectors()):
    assert d1.is_similar_to(d2,
      fast_axis_tolerance=5e-5, slow_axis_tolerance=5e-5, origin_tolerance=5e-5)
  print "OK"
  return

def run():
  if not libtbx.env.has_module("dials_regression"):
    print "Skipping tests in " + __file__ + " as dials_regression not present"
    return

  test1()
  test2()

if __name__ == '__main__':
  from libtbx.utils import show_times_at_exit
  show_times_at_exit()
  run()
