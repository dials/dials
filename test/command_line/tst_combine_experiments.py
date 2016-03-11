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
Test combination of multiple experiments and reflections files.

"""

# python imports
from __future__ import division
import os
import libtbx.load_env # required for libtbx.env.find_in_repositories
from libtbx import easy_run
from libtbx.test_utils import open_tmp_directory
from dxtbx.model.experiment.experiment_list import ExperimentListFactory
from dials.array_family import flex

phil_input = """
  input.experiments={0}/data/sweep_002/experiments.json
  input.reflections={0}/data/sweep_002/reflections.pickle
  input.experiments={0}/data/sweep_003/experiments.json
  input.reflections={0}/data/sweep_003/reflections.pickle
  input.experiments={0}/data/sweep_004/experiments.json
  input.reflections={0}/data/sweep_004/reflections.pickle
  input.experiments={0}/data/sweep_005/experiments.json
  input.reflections={0}/data/sweep_005/reflections.pickle
  input.experiments={0}/data/sweep_006/experiments.json
  input.reflections={0}/data/sweep_006/reflections.pickle
  input.experiments={0}/data/sweep_007/experiments.json
  input.reflections={0}/data/sweep_007/reflections.pickle
  input.experiments={0}/data/sweep_009/experiments.json
  input.reflections={0}/data/sweep_009/reflections.pickle
  input.experiments={0}/data/sweep_011/experiments.json
  input.reflections={0}/data/sweep_011/reflections.pickle
  input.experiments={0}/data/sweep_012/experiments.json
  input.reflections={0}/data/sweep_012/reflections.pickle
  input.experiments={0}/data/sweep_013/experiments.json
  input.reflections={0}/data/sweep_013/reflections.pickle
  input.experiments={0}/data/sweep_014/experiments.json
  input.reflections={0}/data/sweep_014/reflections.pickle
  input.experiments={0}/data/sweep_017/experiments.json
  input.reflections={0}/data/sweep_017/reflections.pickle
  input.experiments={0}/data/sweep_018/experiments.json
  input.reflections={0}/data/sweep_018/reflections.pickle
  input.experiments={0}/data/sweep_019/experiments.json
  input.reflections={0}/data/sweep_019/reflections.pickle
  input.experiments={0}/data/sweep_020/experiments.json
  input.reflections={0}/data/sweep_020/reflections.pickle
  input.experiments={0}/data/sweep_021/experiments.json
  input.reflections={0}/data/sweep_021/reflections.pickle
  input.experiments={0}/data/sweep_022/experiments.json
  input.reflections={0}/data/sweep_022/reflections.pickle
  input.experiments={0}/data/sweep_023/experiments.json
  input.reflections={0}/data/sweep_023/reflections.pickle
  input.experiments={0}/data/sweep_024/experiments.json
  input.reflections={0}/data/sweep_024/reflections.pickle
  input.experiments={0}/data/sweep_025/experiments.json
  input.reflections={0}/data/sweep_025/reflections.pickle
  input.experiments={0}/data/sweep_026/experiments.json
  input.reflections={0}/data/sweep_026/reflections.pickle
  input.experiments={0}/data/sweep_027/experiments.json
  input.reflections={0}/data/sweep_027/reflections.pickle
  input.experiments={0}/data/sweep_028/experiments.json
  input.reflections={0}/data/sweep_028/reflections.pickle
  input.experiments={0}/data/sweep_029/experiments.json
  input.reflections={0}/data/sweep_029/reflections.pickle
  input.experiments={0}/data/sweep_030/experiments.json
  input.reflections={0}/data/sweep_030/reflections.pickle
  input.experiments={0}/data/sweep_031/experiments.json
  input.reflections={0}/data/sweep_031/reflections.pickle
  input.experiments={0}/data/sweep_032/experiments.json
  input.reflections={0}/data/sweep_032/reflections.pickle
  input.experiments={0}/data/sweep_033/experiments.json
  input.reflections={0}/data/sweep_033/reflections.pickle
  input.experiments={0}/data/sweep_035/experiments.json
  input.reflections={0}/data/sweep_035/reflections.pickle
  input.experiments={0}/data/sweep_036/experiments.json
  input.reflections={0}/data/sweep_036/reflections.pickle
  input.experiments={0}/data/sweep_037/experiments.json
  input.reflections={0}/data/sweep_037/reflections.pickle
  input.experiments={0}/data/sweep_038/experiments.json
  input.reflections={0}/data/sweep_038/reflections.pickle
  input.experiments={0}/data/sweep_040/experiments.json
  input.reflections={0}/data/sweep_040/reflections.pickle
  input.experiments={0}/data/sweep_041/experiments.json
  input.reflections={0}/data/sweep_041/reflections.pickle
  input.experiments={0}/data/sweep_042/experiments.json
  input.reflections={0}/data/sweep_042/reflections.pickle
  input.experiments={0}/data/sweep_043/experiments.json
  input.reflections={0}/data/sweep_043/reflections.pickle
  input.experiments={0}/data/sweep_044/experiments.json
  input.reflections={0}/data/sweep_044/reflections.pickle
  input.experiments={0}/data/sweep_046/experiments.json
  input.reflections={0}/data/sweep_046/reflections.pickle
  input.experiments={0}/data/sweep_047/experiments.json
  input.reflections={0}/data/sweep_047/reflections.pickle
  input.experiments={0}/data/sweep_048/experiments.json
  input.reflections={0}/data/sweep_048/reflections.pickle
 """

def test1():

  dials_regression = libtbx.env.find_in_repositories(
    relative_path="dials_regression",
    test=os.path.isdir)

  data_dir = os.path.join(dials_regression, "refinement_test_data",
                          "multi_narrow_wedges")

  # work in a temporary directory
  cwd = os.path.abspath(os.curdir)
  tmp_dir = open_tmp_directory(suffix="tst_combine_experiments_and_reflections")
  os.chdir(tmp_dir)

  input_phil = phil_input.format(data_dir) + """
 reference_from_experiment.beam=0
 reference_from_experiment.scan=0
 reference_from_experiment.goniometer=0
 reference_from_experiment.detector=0
 """

  with open("input.phil","w") as phil_file:
      phil_file.writelines(input_phil)

  cmd = "dials.combine_experiments input.phil"
  #print cmd

  result = easy_run.fully_buffered(command=cmd).raise_if_errors()
  # load results
  exp = ExperimentListFactory.from_json_file("combined_experiments.json",
              check_format=False)
  ref = flex.reflection_table.from_pickle("combined_reflections.pickle")

  # test the experiments
  assert len(exp) == 103
  assert len(exp.crystals()) == 103
  assert len(exp.beams()) == 1
  assert len(exp.scans()) == 1
  assert len(exp.detectors()) == 1
  assert len(exp.goniometers()) == 1
  for e in exp:
    assert e.imageset is not None

  # test the reflections
  assert len(ref) == 11689

  cmd = " ".join([
    "dials.split_experiments",
    "combined_experiments.json",
    "combined_reflections.pickle"])

  result = easy_run.fully_buffered(command=cmd).raise_if_errors()
  for i in range(len(exp)):
    assert os.path.exists("experiments_%03d.json" %i)
    assert os.path.exists("reflections_%03d.pickle" %i)

    exp_single = ExperimentListFactory.from_json_file(
      "experiments_%03d.json" %i, check_format=False)
    ref_single = flex.reflection_table.from_pickle("reflections_%03d.pickle" %i)

    assert len(exp_single) ==1
    assert exp_single[0].crystal == exp[i].crystal
    assert exp_single[0].beam == exp[i].beam
    assert exp_single[0].detector == exp[i].detector
    assert exp_single[0].scan == exp[i].scan
    assert exp_single[0].goniometer == exp[i].goniometer
    assert exp_single[0].imageset == exp[i].imageset
    assert len(ref_single) == len(ref.select(ref['id'] == i))
    assert ref_single['id'].all_eq(0)

  cmd = " ".join([
    "dials.split_experiments",
    "combined_experiments.json",
    "output.experiments_prefix=test"])

  result = easy_run.fully_buffered(command=cmd).raise_if_errors()
  for i in range(len(exp)):
    assert os.path.exists("test_%03d.json" %i)

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
