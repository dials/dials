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
Test refinement of multiple narrow sweeps.

"""

# python imports
from __future__ import division
import os
import shutil
import libtbx.load_env # required for libtbx.env.find_in_repositories
from libtbx import easy_run
from libtbx.test_utils import open_tmp_directory
from scitbx import matrix
from dxtbx.model.experiment.experiment_list import ExperimentListFactory
from dials.array_family import flex

phil_input = """experiments={0}/data/sweep_002/experiments.json \
  reflections={0}/data/sweep_002/reflections.pickle \
  experiments={0}/data/sweep_003/experiments.json \
  reflections={0}/data/sweep_003/reflections.pickle \
  experiments={0}/data/sweep_004/experiments.json \
  reflections={0}/data/sweep_004/reflections.pickle \
  experiments={0}/data/sweep_005/experiments.json \
  reflections={0}/data/sweep_005/reflections.pickle \
  experiments={0}/data/sweep_006/experiments.json \
  reflections={0}/data/sweep_006/reflections.pickle \
  experiments={0}/data/sweep_007/experiments.json \
  reflections={0}/data/sweep_007/reflections.pickle \
  experiments={0}/data/sweep_009/experiments.json \
  reflections={0}/data/sweep_009/reflections.pickle \
  experiments={0}/data/sweep_011/experiments.json \
  reflections={0}/data/sweep_011/reflections.pickle \
  experiments={0}/data/sweep_012/experiments.json \
  reflections={0}/data/sweep_012/reflections.pickle \
  experiments={0}/data/sweep_013/experiments.json \
  reflections={0}/data/sweep_013/reflections.pickle \
  experiments={0}/data/sweep_014/experiments.json \
  reflections={0}/data/sweep_014/reflections.pickle \
  experiments={0}/data/sweep_017/experiments.json \
  reflections={0}/data/sweep_017/reflections.pickle \
  experiments={0}/data/sweep_018/experiments.json \
  reflections={0}/data/sweep_018/reflections.pickle \
  experiments={0}/data/sweep_019/experiments.json \
  reflections={0}/data/sweep_019/reflections.pickle \
  experiments={0}/data/sweep_020/experiments.json \
  reflections={0}/data/sweep_020/reflections.pickle"""

class Test(object):

  def __init__(self):

    dials_regression = libtbx.env.find_in_repositories(
      relative_path="dials_regression",
      test=os.path.isdir)

    self._data_dir = os.path.join(dials_regression, "refinement_test_data",
                            "multi_narrow_wedges")

    # set up a temporary directory
    self._cwd = os.path.abspath(os.curdir)
    self._tmp_dir = open_tmp_directory(suffix="tst_refine_multi_wedges")
    os.chdir(self._tmp_dir)

  def _finish(self):
    """Change back to the original directory and delete the temp dir"""

    os.chdir(self._cwd)
    shutil.rmtree(self._tmp_dir)
    return

  def _combine(self):
    """Combine all the separate sweeps"""

    cmd = "dials.combine_experiments " + phil_input.format(
      self._data_dir) + " reference_from_experiment.beam=0 " + \
      "reference_from_experiment.goniometer=0"+ \
      " reference_from_experiment.detector=0"
    result = easy_run.fully_buffered(command=cmd).raise_if_errors()
    return

  def _refine(self):
    """Do refinement and load the results"""

    # turn off outlier rejection so that test takes about 4s rather than 10s
    cmd = "dials.refine combined_experiments.json combined_reflections.pickle outlier.algorithm=null"
    result = easy_run.fully_buffered(command=cmd).raise_if_errors()

    self._refined_experiments = ExperimentListFactory.from_json_file(
      "refined_experiments.json", check_format=False)
    return

  def run(self):
    """Run the test"""

    try:
      self._combine()
      self._refine()
    finally:
      self._finish()

    return

  def regression(self):
    """Check results are as expected"""

    regression_experiments = ExperimentListFactory.from_json_file(
      os.path.join(self._data_dir, "regression_experiments.json"),
      check_format=False)

    for e1, e2 in zip(self._refined_experiments, regression_experiments):
      assert e1.crystal.is_similar_to(e2.crystal)
      # FIXME need is_similar_to for detector that checks geometry
      #assert e1.detector == e2.detector
      s0_1 = matrix.col(e1.beam.get_unit_s0())
      s0_2 = matrix.col(e1.beam.get_unit_s0())
      assert s0_1.accute_angle(s0_2, deg=True) < 0.0057 # ~0.1 mrad
    print "OK"
    return

def run():
  if not libtbx.env.has_module("dials_regression"):
    print "Skipping tests in " + __file__ + " as dials_regression not present"
    return

  tst = Test()
  tst.run()
  tst.regression()

if __name__ == '__main__':
  from libtbx.utils import show_times_at_exit
  show_times_at_exit()
  run()
