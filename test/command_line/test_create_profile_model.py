from __future__ import absolute_import, division, print_function

import os
import procrunner
import pytest

def test_basic(dials_regression, tmpdir):
  tmpdir.chdir()

  # Call dials.create_profile_model
  result = procrunner.run_process([
      'dials.create_profile_model',
      os.path.join(dials_regression, "integration_test_data", "i04-weak-data2", 'experiments.json'),
      os.path.join(dials_regression, "integration_test_data", "i04-weak-data2", 'indexed.pickle'),
      'sigma_m_algorithm=basic',
  ])
  assert result['exitcode'] == 0
  assert result['stderr'] == ''
  assert os.path.exists('experiments_with_profile_model.json')

  from dxtbx.model.experiment_list import ExperimentListFactory
  experiments =  ExperimentListFactory.from_json_file(
      "experiments_with_profile_model.json",
      check_format=False)
  sigma_b = experiments[0].profile.sigma_b(deg=True)
  sigma_m = experiments[0].profile.sigma_m(deg=True)
  assert sigma_b == pytest.approx(0.02195, abs=1e-3)
  assert sigma_m == pytest.approx(0.06833, abs=1e-3)

def test_extended(dials_regression, tmpdir):
  tmpdir.chdir()

  # Call dials.create_profile_model
  result = procrunner.run_process([
      'dials.create_profile_model',
      os.path.join(dials_regression, "integration_test_data", "i04-weak-data2", 'experiments.json'),
      os.path.join(dials_regression, "integration_test_data", "i04-weak-data2", 'indexed.pickle'),
  ])
  assert result['exitcode'] == 0
  assert result['stderr'] == ''
  assert os.path.exists('experiments_with_profile_model.json')

  from dxtbx.model.experiment_list import ExperimentListFactory
  experiments =  ExperimentListFactory.from_json_file(
      "experiments_with_profile_model.json",
      check_format=False)
  sigma_b = experiments[0].profile.sigma_b(deg=True)
  sigma_m = experiments[0].profile.sigma_m(deg=True)
  assert sigma_b == pytest.approx(0.02195, abs=1e-3)
  assert sigma_m == pytest.approx(0.04187, abs=1e-3)
