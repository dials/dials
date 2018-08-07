"""
Test the command line script dials.scale, for successful completion.
Also tested are the command line files dials.plot_scaling_models and
dials.plot_scaling_outliers
"""

from __future__ import absolute_import, division, print_function

import os
import cPickle as pickle
import pytest
from libtbx import easy_run
from dxtbx.serialize import load

class run_one_scaling(object):
  """Class to run the dials.scale algorithm."""
  def __init__(self, pickle_path_list, sweep_path_list, extra_args, plot=True):
    args = ["dials.scale"] + pickle_path_list + sweep_path_list + extra_args
    command = " ".join(args)
    print(command)
    _ = easy_run.fully_buffered(command=command).raise_if_errors()
    assert os.path.exists("scaled_experiments.json")
    assert os.path.exists("scaled.pickle")

    with open('scaled.pickle', 'rb') as fh:
      table = pickle.load(fh)

    assert 'inverse_scale_factor' in table
    assert 'inverse_scale_factor_variance' in table

    if plot:
      args = ["dials.plot_scaling_models", "scaled.pickle",
        "scaled_experiments.json"]
      command = " ".join(args)
      print(command)
      _ = easy_run.fully_buffered(command=command).raise_if_errors()
      args = ["dials.plot_scaling_outliers", "scaled.pickle",
        "scaled_experiments.json"]
      command = " ".join(args)
      print(command)
      _ = easy_run.fully_buffered(command=command).raise_if_errors()

def get_merging_stats(scaled_unmerged_mtz, anomalous=False, n_bins=20,
                  use_internal_variance=False, eliminate_sys_absent=False,
                  data_labels=None):
  import iotbx.merging_statistics
  i_obs = iotbx.merging_statistics.select_data(
    scaled_unmerged_mtz, data_labels=data_labels)
  i_obs = i_obs.customized_copy(anomalous_flag=False, info=i_obs.info())
  result = iotbx.merging_statistics.dataset_statistics(
    i_obs=i_obs,
    n_bins=n_bins,
    anomalous=anomalous,
    use_internal_variance=use_internal_variance,
    eliminate_sys_absent=eliminate_sys_absent,
  )
  return result

@pytest.mark.dataset_test
def test_scale_physical(dials_regression, tmpdir):
  """Test standard scaling of one dataset."""

  data_dir = os.path.join(dials_regression, "xia2-28",)
  pickle_path = os.path.join(data_dir, "20_integrated.pickle")
  sweep_path = os.path.join(data_dir, "20_integrated_experiments.json")
  extra_args = ["model=physical", "merged_mtz=merged.mtz",
    "unmerged_mtz=unmerged.mtz", "use_free_set=1", "outlier_rejection=simple"]

  with tmpdir.as_cwd():
    _ = run_one_scaling([pickle_path], [sweep_path], extra_args)
    assert os.path.exists("unmerged.mtz")
    assert os.path.exists("merged.mtz")

    # Now inspect output, check it hasn't changed drastically, or if so verify
    # that the new behaviour is more correct and update test accordingly.
    result = get_merging_stats("unmerged.mtz")
    assert result.overall.r_pim < 0.024 #at 07/08/18, value was 0.0234517
    assert result.overall.cc_one_half > 0.9955 # at 07/08/18, value was 0.99597
    assert result.overall.n_obs > 2300 # at 07/08/18, was 2309

  # run again with the concurrent scaling option turned off and the 'standard'
  # outlier rejection
  extra_args = ["model=physical", "merged_mtz=merged.mtz",
    "unmerged_mtz=unmerged.mtz", "use_free_set=1", "outlier_rejection=standard",
    "concurrent=False"]
  with tmpdir.as_cwd():
    _ = run_one_scaling([pickle_path], [sweep_path], extra_args)
    assert os.path.exists("scale_model.png")
    assert os.path.exists("absorption_surface.png")
    assert os.path.exists("outliers.png")

    # Now inspect output, check it hasn't changed drastically, or if so verify
    # that the new behaviour is more correct and update test accordingly.
    result = get_merging_stats("unmerged.mtz")
    assert result.overall.r_pim < 0.035 #at 07/08/18, value was 0.034104
    assert result.overall.cc_one_half > 0.9935 # at 07/08/18, value was 0.99388
    assert result.overall.n_obs > 2310 # at 07/08/18, was 2319

@pytest.mark.dataset_test
def test_scale_optimise_errors(dials_regression, tmpdir):
  """Test standard scaling of one dataset with error optimisation."""
  data_dir = os.path.join(dials_regression, "xia2-28",)
  pickle_path = os.path.join(data_dir, "20_integrated.pickle")
  sweep_path = os.path.join(data_dir, "20_integrated_experiments.json")
  extra_args = ["model=physical", "optimise_errors=True"]
  with tmpdir.as_cwd():
    _ = run_one_scaling([pickle_path], [sweep_path], extra_args)

@pytest.mark.dataset_test
def test_scale_array(dials_regression, tmpdir):
  """Test a standard dataset - ideally needs a large dataset or full matrix
  round may fail. Currently turning off absorption term to avoid
  overparameterisation and failure of full matrix minimisation."""

  data_dir = os.path.join(dials_regression, "xia2-28",)
  pickle_path = os.path.join(data_dir, "20_integrated.pickle")
  sweep_path = os.path.join(data_dir, "20_integrated_experiments.json")
  extra_args = ["model=array", "absorption_term=0", "full_matrix=0"]

  with tmpdir.as_cwd():
    _ = run_one_scaling([pickle_path], [sweep_path], extra_args)
    assert os.path.exists("decay_correction.png")
    assert os.path.exists("outliers.png")

@pytest.mark.dataset_test
def test_multi_scale(dials_regression, tmpdir):
  """Test standard scaling of two datasets."""

  data_dir = os.path.join(dials_regression, "xia2-28",)
  pickle_path_1 = os.path.join(data_dir, "20_integrated.pickle")
  sweep_path_1 = os.path.join(data_dir, "20_integrated_experiments.json")
  pickle_path_2 = os.path.join(data_dir, "25_integrated.pickle")
  sweep_path_2 = os.path.join(data_dir, "25_integrated_experiments.json")
  extra_args = ["unmerged_mtz=unmerged.mtz"]

  with tmpdir.as_cwd():
    _ = run_one_scaling([pickle_path_1, pickle_path_2],
      [sweep_path_1, sweep_path_2], extra_args)
    assert os.path.exists("scale_model_1.png")
    assert os.path.exists("scale_model_2.png")
    assert os.path.exists("absorption_surface_1.png")
    assert os.path.exists("absorption_surface_2.png")
    assert os.path.exists("outliers_1.png")
    assert os.path.exists("outliers_2.png")

    # Now inspect output, check it hasn't changed drastically, or if so verify
    # that the new behaviour is more correct and update test accordingly.
    result = get_merging_stats("unmerged.mtz")
    expected_nobs = 5460
    assert abs(result.overall.n_obs - expected_nobs) < 10
    assert result.overall.r_pim < 0.022 #at 07/08/18, value was 0.021837
    assert result.overall.cc_one_half > 0.9975 # at 07/08/18, value was 0.99810

    #run again, optimising errors, and continuing from where last run left off.
    extra_args = ["optimise_errors=True", "unmerged_mtz=unmerged.mtz"]
    _ = run_one_scaling(["scaled.pickle"], ["scaled_experiments.json"],
      extra_args)
    # Now inspect output, check it hasn't changed drastically, or if so verify
    # that the new behaviour is more correct and update test accordingly.
    # Note: error optimisation currently appears to give worse results here!
    result = get_merging_stats("unmerged.mtz")
    expected_nobs = 5459
    assert abs(result.overall.n_obs - expected_nobs) < 10
    assert result.overall.r_pim < 0.023 #at 07/08/18, value was 0.022722
    assert result.overall.cc_one_half > 0.9965 # at 07/08/18, value was 0.996925

@pytest.mark.dataset_test
def test_targeted_scaling(dials_regression, tmpdir):
  """Test the targeted scaling workflow."""
  data_dir = os.path.join(dials_regression, "xia2-28",)
  pickle_path_1 = os.path.join(data_dir, "20_integrated.pickle")
  sweep_path_1 = os.path.join(data_dir, "20_integrated_experiments.json")
  pickle_path_2 = os.path.join(data_dir, "25_integrated.pickle")
  sweep_path_2 = os.path.join(data_dir, "25_integrated_experiments.json")

  extra_args = ["model=physical"]

  with tmpdir.as_cwd():

    args = ["dials.scale"] + [pickle_path_1] + [sweep_path_1] + extra_args
    command = " ".join(args)
    print(command)
    _ = easy_run.fully_buffered(command=command).raise_if_errors()
    assert os.path.exists("scaled_experiments.json")
    assert os.path.exists("scaled.pickle")

    experiments_list = load.experiment_list(
      "scaled_experiments.json", check_format=False)
    assert len(experiments_list.scaling_models()) == 1

    # Once individual has run, do targeted scaling of second dataset.
    # Use this as a chance to test the KB model as well.
    extra_args = ["model=KB"]
    args = ["dials.scale"] + ["scaled.pickle"] + ["scaled_experiments.json"] +\
      [pickle_path_2] + [sweep_path_2] + extra_args
    command = " ".join(args)
    print(command)
    _ = easy_run.fully_buffered(command=command).raise_if_errors()
    assert os.path.exists("scaled_experiments.json")
    assert os.path.exists("scaled.pickle")

    experiments_list = load.experiment_list(
      "scaled_experiments.json", check_format=False)
    assert len(experiments_list.scaling_models()) == 2
    assert experiments_list.scaling_models()[0].id_ == 'physical'
    assert experiments_list.scaling_models()[1].id_ == 'KB'

    extra_args = ["model=KB", "only_target=True"]
    args = ["dials.scale"] + ["scaled.pickle"] + ["scaled_experiments.json"] +\
      [pickle_path_2] + [sweep_path_2] + extra_args
    command = " ".join(args)
    print(command)
    _ = easy_run.fully_buffered(command=command).raise_if_errors()
    assert os.path.exists("scaled_experiments.json")
    assert os.path.exists("scaled.pickle")
