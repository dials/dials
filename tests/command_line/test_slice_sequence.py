from __future__ import annotations

from pathlib import Path

import procrunner
import pytest

from dxtbx.serialize import load

from dials.array_family import flex


def test_slice_sequence_and_compare_with_expected_results(dials_regression, tmp_path):
    # use the i04_weak_data for this test
    data_dir = Path(dials_regression, "refinement_test_data", "i04_weak_data")
    experiments_path = data_dir / "experiments.json"
    pickle_path = data_dir / "indexed_strong.pickle"

    assert experiments_path.is_file()
    assert pickle_path.is_file()

    result = procrunner.run(
        ["dials.slice_sequence", experiments_path, pickle_path, "image_range=1 20"],
        working_directory=tmp_path,
    )
    assert not result.returncode and not result.stderr

    # load results
    sliced_exp = load.experiment_list(
        tmp_path / "experiments_1_20.expt", check_format=False
    )[0]
    sliced_refs = flex.reflection_table.from_file(tmp_path / "indexed_strong_1_20.refl")

    # simple test of results
    assert sliced_exp.scan.get_image_range() == (1, 20)
    assert len(sliced_refs) == 3670


def test_slice_sequence_with_first_images_missing(dials_regression, tmp_path):
    """Test slicing where scan image range does not start at 1, exercising
    a case that exposed a bug"""

    # use the i04_weak_data for this test
    data_dir = Path(dials_regression, "refinement_test_data", "i04_weak_data")
    experiments_path = data_dir / "experiments.json"

    # first slice
    result = procrunner.run(
        ["dials.slice_sequence", experiments_path, "image_range=5,20"],
        working_directory=tmp_path,
    )
    assert not result.returncode and not result.stderr

    # second slice
    result = procrunner.run(
        ["dials.slice_sequence", "experiments_5_20.expt", "image_range=10,20"],
        working_directory=tmp_path,
    )
    assert not result.returncode and not result.stderr

    sliced_exp = load.experiment_list(
        tmp_path / "experiments_5_20_10_20.expt", check_format=False
    )[0]
    assert sliced_exp.scan.get_image_range() == (10, 20)
    assert sliced_exp.scan.get_array_range() == (9, 20)
    assert sliced_exp.scan.get_oscillation()[0] == pytest.approx(83.35)


def test_slice_sequence_to_degree_blocks(dials_data, tmp_path):
    """Slice data into 10 degree blocks i.e. 17 datasets"""
    expt = dials_data("l_cysteine_4_sweeps_scaled", pathlib=True) / "scaled_30.expt"
    refl = dials_data("l_cysteine_4_sweeps_scaled", pathlib=True) / "scaled_30.refl"
    procrunner.run(
        [
            "dials.slice_sequence",
            "block_size=10",
            "output.experiments=sliced.expt",
            "output.reflections=sliced.refl",
            expt,
            refl,
        ],
        working_directory=tmp_path,
    )

    sliced_expts = load.experiment_list(tmp_path / "sliced.expt", check_format=False)
    assert len(sliced_expts) == 17
    sliced_refl = flex.reflection_table.from_file(tmp_path / "sliced.refl")
    assert len(set(sliced_refl.experiment_identifiers().values())) == 17
    sliced_refl.assert_experiment_identifiers_are_consistent(sliced_expts)


def test_slice_sequence_with_scan_varying_crystal(dials_data, tmp_path):
    """test slicing keeps a scan-varying crystal"""

    expt = dials_data("l_cysteine_4_sweeps_scaled", pathlib=True) / "scaled_30.expt"
    procrunner.run(
        [
            "dials.slice_sequence",
            "image_range=10,20",
            "output.experiments=sliced.expt",
            expt,
        ],
        working_directory=tmp_path,
    )
    orig = load.experiment_list(expt, check_format=False)[0]

    sliced = load.experiment_list(tmp_path / "sliced.expt", check_format=False)[0]

    assert sliced.crystal.num_scan_points == 12

    orig_UB = [
        orig.crystal.get_A_at_scan_point(i) for i in range(orig.crystal.num_scan_points)
    ]
    sliced_UB = [
        sliced.crystal.get_A_at_scan_point(i)
        for i in range(sliced.crystal.num_scan_points)
    ]

    for a, b in zip(orig_UB[9:21], sliced_UB):
        assert a == pytest.approx(b)
