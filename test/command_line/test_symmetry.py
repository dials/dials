from __future__ import absolute_import, division, print_function

import os
import pytest
import procrunner


def test_symmetry(dials_regression, run_in_tmpdir):
  """Simple test to check that dials.symmetry completes"""

  result = procrunner.run_process([
    'dials.symmetry',
    os.path.join(dials_regression, "xia2-28", "20_integrated_experiments.json"),
    os.path.join(dials_regression, "xia2-28", "20_integrated.pickle"),
    os.path.join(dials_regression, "xia2-28", "25_integrated_experiments.json"),
    os.path.join(dials_regression, "xia2-28", "25_integrated.pickle")
  ])
  assert result['exitcode'] == 0
  assert result['stderr'] == ''
  assert os.path.exists("reindexed_reflections.pickle")
  assert os.path.exists("reindexed_experiments.json")
