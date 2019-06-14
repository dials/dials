from __future__ import absolute_import, division, print_function

import six.moves.cPickle as pickle
import pytest
import os
from libtbx import easy_run
from dxtbx.model.experiment_list import ExperimentListFactory


def test_slice_sweep_and_compare_with_expected_results(dials_regression, run_in_tmpdir):
    # use the i04_weak_data for this test
    data_dir = os.path.join(dials_regression, "refinement_test_data", "i04_weak_data")
    experiments_path = os.path.join(data_dir, "experiments.json")
    pickle_path = os.path.join(data_dir, "indexed_strong.pickle")

    for pth in (experiments_path, pickle_path):
        assert os.path.exists(pth)

    cmd = (
        "dials.slice_sweep "
        + experiments_path
        + " "
        + pickle_path
        + ' "image_range=1 20"'
    )

    result = easy_run.fully_buffered(command=cmd).raise_if_errors()
    # load results
    sliced_exp = ExperimentListFactory.from_json_file(
        "experiments_1_20.expt", check_format=False
    )[0]
    with open("indexed_strong_1_20.refl", "rb") as f:
        sliced_refs = pickle.load(f)

    # simple test of results
    assert sliced_exp.scan.get_image_range() == (1, 20)
    assert len(sliced_refs) == 3670


def test_slice_sweep_with_first_images_missing(dials_regression, run_in_tmpdir):
    """Test slicing where scan image range does not start at 1, exercising
    a case that exposed a bug"""

    # use the i04_weak_data for this test
    data_dir = os.path.join(dials_regression, "refinement_test_data", "i04_weak_data")
    experiments_path = os.path.join(data_dir, "experiments.json")

    # first slice
    cmd = "dials.slice_sweep " + experiments_path + " image_range=5,20"
    result = easy_run.fully_buffered(command=cmd).raise_if_errors()

    # second slice
    cmd = "dials.slice_sweep experiments_5_20.expt image_range=10,20"
    result = easy_run.fully_buffered(command=cmd).raise_if_errors()

    sliced_exp = ExperimentListFactory.from_json_file(
        "experiments_5_20_10_20.expt", check_format=False
    )[0]
    assert sliced_exp.scan.get_image_range() == (10, 20)
    assert sliced_exp.scan.get_array_range() == (9, 20)
    assert sliced_exp.scan.get_oscillation()[0] == pytest.approx(83.35)
