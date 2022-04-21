from __future__ import annotations

import glob
import os

import pytest

import scitbx
from cctbx import uctbx
from dxtbx.model import ExperimentList
from dxtbx.serialize import load

from dials.command_line import search_beam_position

from ..algorithms.indexing.test_index import run_indexing


def test_search_i04_weak_data_image_range(mocker, run_in_tmp_path, dials_regression):
    """Perform a beam-centre search and check that the output is sane."""

    data_dir = os.path.join(dials_regression, "indexing_test_data", "i04_weak_data")
    reflection_file = os.path.join(data_dir, "full.pickle")
    experiments_file = os.path.join(data_dir, "experiments_import.json")

    args = [
        experiments_file,
        reflection_file,
        "image_range=1,10",
        "image_range=251,260",
        "image_range=531,540",
        "n_macro_cycles=4",
    ]
    from rstbx.indexing_api import dps_extended

    mocker.spy(dps_extended, "get_new_detector")
    search_beam_position.run(args)
    # Check that the last call to get_new_detector was with an offset of close to zero.
    # The final call was to apply the "best" shift to the detector model before
    # returning the updated experiments.
    assert dps_extended.get_new_detector.call_args[0][1].elems == pytest.approx(
        (0, 0, 0), abs=3e-2
    )
    assert os.path.exists("optimised.expt")

    # Compare the shifts between the start and final detector models
    experiments = load.experiment_list(experiments_file, check_format=False)
    optimised_experiments = load.experiment_list("optimised.expt", check_format=False)
    detector_1 = experiments[0].detector
    detector_2 = optimised_experiments[0].detector
    shift = scitbx.matrix.col(detector_1[0].get_origin()) - scitbx.matrix.col(
        detector_2[0].get_origin()
    )
    assert shift.elems == pytest.approx((0.27, -0.12, 0.0), abs=1e-1)


def test_search_multiple(run_in_tmp_path, dials_regression):
    """Perform a beam-centre search and check that the output is sane.

    Do the following:
    1. Run dials.search_beam_centre on two datablocks and two pickled
    reflection tables, as output by dials.find_spots;
      a) Check that the program exits correctly;
      b) Check that it produces the expected output datablock.
    2. Check that the beam centre search has resulted in the expected shift
    in detector origin.
    """

    data_dir = os.path.join(dials_regression, "indexing_test_data", "trypsin")
    pickle_path1 = os.path.join(data_dir, "strong_P1_X6_1_0-1.pickle")
    pickle_path2 = os.path.join(data_dir, "strong_P1_X6_2_0-1.pickle")
    experiments_path1 = os.path.join(data_dir, "datablock_P1_X6_1.json")
    experiments_path2 = os.path.join(data_dir, "datablock_P1_X6_2.json")

    args = [experiments_path1, experiments_path2, pickle_path1, pickle_path2]
    search_beam_position.run(args)
    assert os.path.exists("optimised.expt")

    experiments = load.experiment_list(experiments_path1, check_format=False)
    optimised_experiments = load.experiment_list("optimised.expt", check_format=False)
    detector_1 = experiments[0].detector
    detector_2 = optimised_experiments[0].detector
    shift = scitbx.matrix.col(detector_1[0].get_origin()) - scitbx.matrix.col(
        detector_2[0].get_origin()
    )
    assert shift.elems == pytest.approx((-0.518, 0.192, 0.0), abs=1e-1)


def test_index_after_search(dials_data, run_in_tmp_path):
    """Integrate the beam centre search with the rest of the toolchain

    Do the following:
    1. Take a known good experiment and perturbate the beam centre
    2. Run dials.search_beam_centre on the perturbated beam centre and original
    reflection table, check for expected output;
    3. Run dials.index with the found beam centre and check that the expected
    unit cell is obtained and that the RMSDs are smaller than or equal to some
    expected values."""

    insulin = dials_data("insulin_processed", pathlib=True)

    # load the original experiment and perturbate the beam centre by a small offset
    experiments = load.experiment_list(insulin / "imported.expt", check_format=False)
    original_origin = experiments[0].detector.hierarchy().get_origin()
    shifted_origin = (
        original_origin[0] - 1.3,
        original_origin[1] + 1.5,
        original_origin[2],
    )
    experiments[0].detector.hierarchy().set_local_frame(
        experiments[0].detector.hierarchy().get_fast_axis(),
        experiments[0].detector.hierarchy().get_slow_axis(),
        shifted_origin,
    )
    assert experiments[0].detector.hierarchy().get_origin() == shifted_origin
    experiments.as_file(run_in_tmp_path / "shifted.expt")

    # search the beam centre
    search_beam_position.run(
        [
            str(run_in_tmp_path / "shifted.expt"),
            str(insulin / "strong.refl"),
        ]
    )
    assert run_in_tmp_path.joinpath("optimised.expt").is_file()

    # check we can actually index the resulting optimized experiments
    expected_unit_cell = uctbx.unit_cell(
        (67.655, 67.622, 67.631, 109.4583, 109.4797, 109.485)
    )
    expected_rmsds = (0.3, 0.3, 0.005)
    expected_hall_symbol = " P 1"
    run_indexing(
        insulin / "strong.refl",
        run_in_tmp_path / "optimised.expt",
        run_in_tmp_path,
        [],
        expected_unit_cell,
        expected_rmsds,
        expected_hall_symbol,
    )


def test_search_single(run_in_tmp_path, dials_regression):
    """Perform a beam-centre search and check that the output is sane.

    Do the following:
    1. Run dials.search_beam_centre on a single datablock and pickled
    reflection table, as output by dials.find_spots;
      a) Check that the program exits correctly;
      b) Check that it produces the expected output datablock.
    2. Check that the beam centre search has resulted in the expected shift
    in detector origin.
    """

    data_dir = os.path.join(dials_regression, "indexing_test_data", "phi_scan")
    pickle_path = os.path.join(data_dir, "strong.pickle")
    experiments_path = os.path.join(data_dir, "datablock.json")

    search_beam_position.run([experiments_path, pickle_path])
    assert run_in_tmp_path.joinpath("optimised.expt").is_file()

    experiments = load.experiment_list(experiments_path, check_format=False)
    original_imageset = experiments.imagesets()[0]
    optimized_experiments = load.experiment_list("optimised.expt", check_format=False)
    detector_1 = original_imageset.get_detector()
    detector_2 = optimized_experiments.detectors()[0]
    shift = scitbx.matrix.col(detector_1[0].get_origin()) - scitbx.matrix.col(
        detector_2[0].get_origin()
    )
    assert shift.elems == pytest.approx((-0.976, 2.497, 0.0), abs=1e-1)


def test_search_small_molecule(dials_data, run_in_tmp_path):
    """Perform a beam-centre search on a multi-sequence data set..

    Do the following:
    1. Run dials.search_beam_centre on a single datablock and pickled
    reflection table containing multiple experiment IDs, as output by
    dials.find_spots;
      a) Check that the program exits correctly;
      b) Check that it produces the expected output datablock.
    2. Check that the beam centre search has resulted in the expected shift
    in detector origin.
    """

    data = dials_data("l_cysteine_dials_output", pathlib=True)
    experiments_path = data / "imported.expt"
    refl_path = data / "strong.refl"

    search_beam_position.run([os.fspath(experiments_path), os.fspath(refl_path)])
    assert run_in_tmp_path.joinpath("optimised.expt").is_file()

    experiments = load.experiment_list(experiments_path, check_format=False)
    optimised_experiments = load.experiment_list("optimised.expt", check_format=False)
    for old_expt, new_expt in zip(experiments, optimised_experiments):
        # assert that the detector fast/slow axes are unchanged from the input experiments
        # the last experiment actually does have a different detector model
        assert (
            old_expt.detector[0].get_slow_axis() == new_expt.detector[0].get_slow_axis()
        )
        assert (
            old_expt.detector[0].get_fast_axis() == new_expt.detector[0].get_fast_axis()
        )
        shift = scitbx.matrix.col(
            old_expt.detector[0].get_origin()
        ) - scitbx.matrix.col(new_expt.detector[0].get_origin())
        assert shift.elems == pytest.approx((0.091, -1.11, 0), abs=1e-2)


def test_multi_sweep_fixed_rotation(dials_regression, run_in_tmp_path):
    data_dir = os.path.join(dials_regression, "indexing_test_data", "multi_sweep")
    reflection_files = sorted(
        glob.glob(os.path.join(data_dir, "SWEEP[1,2]", "index", "*_strong.pickle"))
    )
    experiment_files = sorted(
        glob.glob(
            os.path.join(data_dir, "SWEEP[1,2]", "index", "*_datablock_import.json")
        )
    )

    search_beam_position.run(reflection_files + experiment_files)
    assert run_in_tmp_path.joinpath("optimised.expt").is_file()

    experiments = ExperimentList()
    for path in experiment_files:
        experiments.extend(load.experiment_list(path, check_format=False))

    optimised_experiments = load.experiment_list("optimised.expt", check_format=False)
    for orig_expt, new_expt in zip(experiments, optimised_experiments):
        shift = scitbx.matrix.col(
            orig_expt.detector[0].get_origin()
        ) - scitbx.matrix.col(new_expt.detector[0].get_origin())
        print(shift)
        assert shift.elems == pytest.approx((2.293, -0.399, 0), abs=1e-2)
