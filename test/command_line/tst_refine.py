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

  # use the i04_weak_data for this test
  data_dir = os.path.join(dials_regression, "refinement_test_data", "i04_weak_data")
  experiments_path = os.path.join(data_dir, "experiments.json")
  pickle_path = os.path.join(data_dir, "indexed_strong.pickle")

  for pth in (experiments_path, pickle_path):
    assert os.path.exists(pth)

  # set close_to_spindle_cutoff to old default
  cmd = "dials.refine close_to_spindle_cutoff=0.05 use_all_reflections=false " + \
        experiments_path + " " + pickle_path
  print cmd

  # work in a temporary directory
  cwd = os.path.abspath(os.curdir)
  tmp_dir = open_tmp_directory(suffix="test_dials_refine")
  os.chdir(tmp_dir)
  try:
    result = easy_run.fully_buffered(command=cmd).raise_if_errors()
    # load results
    reg_exp = ExperimentListFactory.from_json_file(
                os.path.join(data_dir, "regression_experiments.json"),
                check_format=False)[0]
    ref_exp = ExperimentListFactory.from_json_file("refined_experiments.json",
                check_format=False)[0]
  finally:
    os.chdir(cwd)
    # clean up tmp dir
    shutil.rmtree(tmp_dir)

  # test refined models against expected
  assert reg_exp.crystal == ref_exp.crystal
  assert reg_exp.detector == ref_exp.detector
  assert reg_exp.beam == ref_exp.beam

  print "OK"
  return

def test2():
  """Run scan-varying refinement, comparing RMSD table with expected values.
  This test automates what was manually done periodically and recorded in
  dials_regression/refinement_test_data/centroid/README.txt"""

  dials_regression = libtbx.env.find_in_repositories(
    relative_path="dials_regression",
    test=os.path.isdir)

  # use the i04_weak_data for this test
  data_dir = os.path.join(dials_regression, "refinement_test_data", "centroid")
  experiments_path = os.path.join(data_dir, "experiments_XPARM_REGULARIZED.json")
  pickle_path = os.path.join(data_dir, "spot_all_xds.pickle")

  for pth in (experiments_path, pickle_path):
    assert os.path.exists(pth)

  # scan-static refinement first to get refined_experiments.json as start point
  cmd1 = "dials.refine " + experiments_path + " " + pickle_path + \
    " reflections_per_degree=50 use_all_reflections=false " + \
    " outlier.algorithm=null close_to_spindle_cutoff=0.05"
  cmd2 = "dials.refine refined_experiments.json " + pickle_path + \
    " scan_varying=true output.history=history.pickle " + \
    " reflections_per_degree=50 use_all_reflections=false " + \
    " outlier.algorithm=null close_to_spindle_cutoff=0.05"

  # work in a temporary directory
  cwd = os.path.abspath(os.curdir)
  tmp_dir = open_tmp_directory(suffix="test_dials_refine")
  os.chdir(tmp_dir)
  try:
    print cmd1
    result1 = easy_run.fully_buffered(command=cmd1).raise_if_errors()
    print cmd2
    result2 = easy_run.fully_buffered(command=cmd2).raise_if_errors()
    # load and check results
    import cPickle as pickle
    history=pickle.load(open("history.pickle", "r"))

    expected_rmsds = [(0.088488398, 0.114583571, 0.001460382),
                      (0.080489334, 0.086406517, 0.001284069),
                      (0.078835086, 0.086052630, 0.001195882),
                      (0.077476911, 0.086194611, 0.001161143),
                      (0.076755840, 0.086090630, 0.001157239),
                      (0.076586376, 0.085939462, 0.001155641),
                      (0.076603722, 0.085878953, 0.001155065),
                      (0.076611382, 0.085862959, 0.001154863),
                      (0.076608732, 0.085856935, 0.001154384),
                      (0.076605731, 0.085852271, 0.001153858),
                      (0.076604576, 0.085852318, 0.001153643),
                      (0.076603981, 0.085854175, 0.001153594)]

    assert approx_equal(history['rmsd'], expected_rmsds)

    # check that the used_in_refinement flag got set correctly
    rt = flex.reflection_table.from_pickle('refined.pickle')
    uir = rt.get_flags(rt.flags.used_in_refinement)
    assert uir.count(True) == history['num_reflections'][-1]

  finally:
    os.chdir(cwd)
    # clean up tmp dir
    shutil.rmtree(tmp_dir)

  print "OK"
  return


def run():
  if not libtbx.env.has_module("dials_regression"):
    print "Skipping tests in " + __file__ + " as dials_regression not present"
    return

  test1()
  test2()

if __name__ == '__main__':
  from dials.test import cd_auto
  with cd_auto(__file__):
    from libtbx.utils import show_times_at_exit
    show_times_at_exit()
    run()
