from __future__ import absolute_import, division, print_function

import os

import procrunner
import pytest
import scitbx

from dxtbx.serialize import load


def test_search_i04_weak_data_image_range(run_in_tmpdir, dials_regression):
    """Perform a beam-centre search and check that the output is sane."""

    data_dir = os.path.join(dials_regression, "indexing_test_data", "i04_weak_data")
    reflection_file = os.path.join(data_dir, "full.pickle")
    experiments_file = os.path.join(data_dir, "experiments_import.json")

    args = [
        "dials.search_beam_position",
        experiments_file,
        reflection_file,
        "image_range=1,10",
        "image_range=251,260",
        "image_range=531,540",
    ]

    print(args)
    result = procrunner.run(args)
    assert not result.returncode and not result.stderr
    assert os.path.exists("optimised.expt")

    experiments = load.experiment_list(experiments_file, check_format=False)
    original_imageset = experiments.imagesets()[0]
    optimized_experiments = load.experiment_list("optimised.expt", check_format=False)
    detector_1 = original_imageset.get_detector()
    detector_2 = optimized_experiments.detectors()[0]
    shift = scitbx.matrix.col(detector_1[0].get_origin()) - scitbx.matrix.col(
        detector_2[0].get_origin()
    )
    assert shift.elems == pytest.approx((0.27, -0.12, 0.0), abs=1e-1)


def test_search_multiple(run_in_tmpdir, dials_regression):
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

    args = [
        "dials.search_beam_position",
        experiments_path1,
        experiments_path2,
        pickle_path1,
        pickle_path2,
    ]

    print(args)
    result = procrunner.run(args)
    assert not result.returncode and not result.stderr
    assert os.path.exists("optimised.expt")

    experiments = load.experiment_list(experiments_path1, check_format=False)
    original_imageset = experiments.imagesets()[0]
    optimized_experiments = load.experiment_list("optimised.expt", check_format=False)
    detector_1 = original_imageset.get_detector()
    detector_2 = optimized_experiments.detectors()[0]
    shift = scitbx.matrix.col(detector_1[0].get_origin()) - scitbx.matrix.col(
        detector_2[0].get_origin()
    )
    assert shift.elems == pytest.approx((0.037, 0.061, 0.0), abs=1e-1)


def test_index_after_search(dials_data, tmpdir):
    """Integrate the beam centre search with the rest of the toolchain

    Do the following:
    1. Run dials.import with a specified beam centre, check for expected output;
    2. Run dials.find_spots, check for expected output;
    3. Run dials.search_beam_centre on the resultant datablock and pickled
    reflection table, check for expected output;
    4. Run dials.index, using the datablock from the beam centre search,
    and check that the expected unit ecll is obtained and that the RMSDs are
    smaller than or equal to some expected values."""

    dials_data = dials_data("thaumatin_i04").listdir(sort=True)
    g = [f.strpath for f in dials_data if f.ext == ".cbf"]

    # beam centre from image headers: 205.28,210.76 mm
    args = ["dials.import", "mosflm_beam_centre=207,212"] + g
    print(args)
    if os.name != "nt":
        result = procrunner.run(args, working_directory=tmpdir)
        assert not result.returncode and not result.stderr
    else:
        # Can't run this command on Windows,
        # as it will exceed the maximum Windows command length limits.
        # So, instead:
        import mock
        import sys

        with tmpdir.as_cwd():
            with mock.patch.object(sys, "argv", args):
                import dials.command_line.dials_import

                dials.command_line.dials_import.Script().run()
    assert tmpdir.join("imported.expt").check()

    # spot-finding, just need a subset of the data
    args = [
        "dials.find_spots",
        "imported.expt",
        "scan_range=1,10",
        "scan_range=531,540",
    ]
    print(args)
    result = procrunner.run(args, working_directory=tmpdir)
    assert not result.returncode and not result.stderr
    assert tmpdir.join("strong.refl").check()

    # actually run the beam centre search
    args = ["dials.search_beam_position", "imported.expt", "strong.refl"]
    print(args)
    result = procrunner.run(args, working_directory=tmpdir)
    assert not result.returncode and not result.stderr
    assert tmpdir.join("optimised.expt").check()

    # look at the results
    from dxtbx.serialize import load

    experiments = load.experiment_list(
        tmpdir.join("imported.expt").strpath, check_format=False
    )
    original_imageset = experiments.imagesets()[0]
    optimized_experiments = load.experiment_list(
        tmpdir.join("optimised.expt").strpath, check_format=False
    )
    detector_1 = original_imageset.get_detector()
    detector_2 = optimized_experiments.detectors()[0]
    shift = scitbx.matrix.col(detector_1[0].get_origin()) - scitbx.matrix.col(
        detector_2[0].get_origin()
    )
    print(shift)

    # check we can actually index the resulting optimized experiments
    from cctbx import uctbx
    from dials.algorithms.indexing.test_index import run_indexing

    expected_unit_cell = uctbx.unit_cell(
        (57.780, 57.800, 150.017, 89.991, 89.990, 90.007)
    )
    expected_rmsds = (0.06, 0.05, 0.001)
    expected_hall_symbol = " P 1"
    run_indexing(
        "strong.refl",
        "optimised.expt",
        tmpdir,
        [],
        expected_unit_cell,
        expected_rmsds,
        expected_hall_symbol,
    )


def test_search_single(run_in_tmpdir, dials_regression):
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

    args = ["dials.search_beam_position", experiments_path, pickle_path]
    print(args)
    result = procrunner.run(args)
    assert not result.returncode and not result.stderr
    assert os.path.exists("optimised.expt")

    from dxtbx.serialize import load

    experiments = load.experiment_list(experiments_path, check_format=False)
    original_imageset = experiments.imagesets()[0]
    optimized_experiments = load.experiment_list("optimised.expt", check_format=False)
    detector_1 = original_imageset.get_detector()
    detector_2 = optimized_experiments.detectors()[0]
    shift = scitbx.matrix.col(detector_1[0].get_origin()) - scitbx.matrix.col(
        detector_2[0].get_origin()
    )
    assert shift.elems == pytest.approx((-0.976, 2.497, 0.0), abs=1e-1)


def test_search_small_molecule(dials_data, run_in_tmpdir):
    """Perform a beam-centre search on a multi-sweep data set..

    Do the following:
    1. Run dials.search_beam_centre on a single datablock and pickled
    reflection table containing multiple experiment IDs, as output by
    dials.find_spots;
      a) Check that the program exits correctly;
      b) Check that it produces the expected output datablock.
    2. Check that the beam centre search has resulted in the expected shift
    in detector origin.
    """

    data = dials_data("l_cysteine_dials_output")
    datablock_path = data.join("datablock.json").strpath
    pickle_path = data.join("strong.pickle").strpath

    args = ["dials.search_beam_position", datablock_path, pickle_path]
    print(args)
    result = procrunner.run(args)
    assert not result.returncode and not result.stderr
    assert os.path.exists("optimised.expt")

    from dxtbx.serialize import load

    datablocks = load.datablock(datablock_path, check_format=False)
    original_imageset = datablocks[0].extract_imagesets()[0]
    detector_1 = original_imageset.get_detector()
    optimized_experiments = load.experiment_list("optimised.expt", check_format=False)
    detector_2 = optimized_experiments[0].detector
    shift = scitbx.matrix.col(detector_1[0].get_origin()) - scitbx.matrix.col(
        detector_2[0].get_origin()
    )
    print(shift)
    assert shift.elems == pytest.approx((0.11, -1.03, 0.0), abs=1e-1)
