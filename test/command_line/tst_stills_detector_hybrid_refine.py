#!/usr/bin/env cctbx.python

#
#  Copyright (C) (2015) STFC Rutherford Appleton Laboratory, UK.
#
#  Author: David Waterman.
#
#  This code is distributed under the BSD license, a copy of which is
#  included in the root directory of this package.
#

"""
Test dials.stills_detector_hybrid_refine by running a short job

"""

# python imports
from __future__ import division
import os
import shutil
import libtbx.load_env # required for libtbx.env.find_in_repositories
from libtbx import easy_run
from libtbx.test_utils import open_tmp_directory

def test1(averaged_reference_detector=False):

  dials_regression = libtbx.env.find_in_repositories(
    relative_path="dials_regression",
    test=os.path.isdir)

  # use 20 indexed pickles from CXI for this test
  data_dir = os.path.join(dials_regression, "stills_test_data",
    "cspad_indexing_results", "cxid9114_r0097")

  experiments_paths = ["idx-20140615231407812_refined_experiments.json",
                       "idx-20140615231408153_refined_experiments.json",
                       "idx-20140615231408912_refined_experiments.json",
                       "idx-20140615231409020_refined_experiments.json",
                       "idx-20140615231409470_refined_experiments.json",
                       "idx-20140615231410120_refined_experiments.json",
                       "idx-20140615231411153_refined_experiments.json",
                       "idx-20140615231411170_refined_experiments.json",
                       "idx-20140615231411970_refined_experiments.json",
                       "idx-20140615231412103_refined_experiments.json",
                       "idx-20140615231412495_refined_experiments.json",
                       "idx-20140615231413370_refined_experiments.json",
                       "idx-20140615231413878_refined_experiments.json",
                       "idx-20140615231414128_refined_experiments.json",
                       "idx-20140615231414461_refined_experiments.json",
                       "idx-20140615231415353_refined_experiments.json",
                       "idx-20140615231415761_refined_experiments.json",
                       "idx-20140615231415819_refined_experiments.json",
                       "idx-20140615231416328_refined_experiments.json",
                       "idx-20140615231416694_refined_experiments.json"]

  reflections_paths = ["idx-20140615231407812_indexed.pickle",
                       "idx-20140615231408153_indexed.pickle",
                       "idx-20140615231408912_indexed.pickle",
                       "idx-20140615231409020_indexed.pickle",
                       "idx-20140615231409470_indexed.pickle",
                       "idx-20140615231410120_indexed.pickle",
                       "idx-20140615231411153_indexed.pickle",
                       "idx-20140615231411170_indexed.pickle",
                       "idx-20140615231411970_indexed.pickle",
                       "idx-20140615231412103_indexed.pickle",
                       "idx-20140615231412495_indexed.pickle",
                       "idx-20140615231413370_indexed.pickle",
                       "idx-20140615231413878_indexed.pickle",
                       "idx-20140615231414128_indexed.pickle",
                       "idx-20140615231414461_indexed.pickle",
                       "idx-20140615231415353_indexed.pickle",
                       "idx-20140615231415761_indexed.pickle",
                       "idx-20140615231415819_indexed.pickle",
                       "idx-20140615231416328_indexed.pickle",
                       "idx-20140615231416694_indexed.pickle"]

  cmd = "dev.dials.stills_detector_hybrid_refine "
  for exp in experiments_paths:
    exp = os.path.join(data_dir, exp)
    cmd += "experiments={0} ".format(exp)
  for ref in reflections_paths:
    ref = os.path.join(data_dir, ref)
    cmd += "reflections={0} ".format(ref)

  if averaged_reference_detector:
    # specify hierarchy_level=0
    cmd += "detector_phase.refinement.parameterisation.detector.hierarchy_level=0 "
    cmd += "reference_detector=average"
  else:
    # specify hierarchy_level=1
    cmd += "detector_phase.refinement.parameterisation.detector.hierarchy_level=1"

  # work in a temporary directory
  cwd = os.path.abspath(os.curdir)
  tmp_dir = open_tmp_directory(suffix="tst_stills_detector_hybrid_refine")
  os.chdir(tmp_dir)
  try:
    result = easy_run.fully_buffered(command=cmd).raise_if_errors()
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
  test1(averaged_reference_detector=True)

if __name__ == '__main__':
  from dials.test import cd_auto
  with cd_auto(__file__):
    from libtbx.utils import show_times_at_exit
    show_times_at_exit()
    run()
