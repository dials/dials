"""
Test refinement of multiple narrow sweeps.
"""

from __future__ import absolute_import, division, print_function

import os

import procrunner
from libtbx.test_utils import open_tmp_directory
from scitbx import matrix
from dxtbx.model.experiment_list import ExperimentListFactory

def test(dials_regression, tmpdir):
  tmpdir.chdir()

  data_dir = os.path.join(dials_regression, "refinement_test_data", "multi_narrow_wedges")

  selection = (2,3,4,5,6,7,9,11,12,13,14,17,18,19,20)

  # Combine all the separate sweeps

  result = procrunner.run_process([
      "dials.combine_experiments",
      "reference_from_experiment.beam=0",
      "reference_from_experiment.goniometer=0",
      "reference_from_experiment.detector=0",
  ] + [
      "experiments={0}/data/sweep_%03d/experiments.json".format(data_dir) % n for n in selection
  ] + [
      "reflections={0}/data/sweep_%03d/reflections.pickle".format(data_dir) % n for n in selection
  ])
  assert result['exitcode'] == 0
  assert result['stderr'] == ''

  # Do refinement and load the results

  # turn off outlier rejection so that test takes about 4s rather than 10s
  # set close_to_spindle_cutoff to old default
  result = procrunner.run_process([
      "dials.refine",
      "combined_experiments.json",
      "combined_reflections.pickle",
      "outlier.algorithm=null",
      "close_to_spindle_cutoff=0.05",
  ])
  assert result['exitcode'] == 0
  assert result['stderr'] == ''

  refined_experiments = ExperimentListFactory.from_json_file("refined_experiments.json", check_format=False)

  # Check results are as expected

  regression_experiments = ExperimentListFactory.from_json_file(
      os.path.join(data_dir, "regression_experiments.json"),
      check_format=False)

  for e1, e2 in zip(refined_experiments, regression_experiments):
    assert e1.crystal.is_similar_to(e2.crystal)
    # FIXME need is_similar_to for detector that checks geometry
    #assert e1.detector == e2.detector
    s0_1 = matrix.col(e1.beam.get_unit_s0())
    s0_2 = matrix.col(e1.beam.get_unit_s0())
    assert s0_1.accute_angle(s0_2, deg=True) < 0.0057 # ~0.1 mrad
