from __future__ import absolute_import, division, print_function

import json
import os

import procrunner
import pytest
from cctbx import uctbx
from dxtbx.serialize import load


def test_refine_bravais_settings(dials_regression, tmpdir):
    data_dir = os.path.join(dials_regression, "indexing_test_data", "i04_weak_data")
    pickle_path = os.path.join(data_dir, "indexed.pickle")
    experiments_path = os.path.join(data_dir, "experiments.json")
    result = procrunner.run(
        [
            "dials.refine_bravais_settings",
            pickle_path,
            experiments_path,
            "reflections_per_degree=5",
            "minimum_sample_size=500",
            "beam.fix=all",
            "detector.fix=all",
            "prefix=tst_",
        ],
        working_directory=tmpdir,
    )
    assert not result.returncode and not result.stderr
    for i in range(1, 10):
        assert tmpdir.join("tst_bravais_setting_%i.expt" % i).check()

    experiments_list = load.experiment_list(
        tmpdir.join("tst_bravais_setting_9.expt").strpath, check_format=False
    )
    assert len(experiments_list) == 1
    assert (
        experiments_list[0]
        .crystal.get_unit_cell()
        .is_similar_to(uctbx.unit_cell((57.782, 57.782, 150.011, 90, 90, 90)))
    )
    assert experiments_list[0].crystal.get_space_group().type().hall_symbol() == " P 4"

    assert tmpdir.join("tst_bravais_summary.json").check()
    with tmpdir.join("tst_bravais_summary.json").open("rb") as fh:
        bravais_summary = json.load(fh)
    assert set(bravais_summary) == {"1", "2", "3", "4", "5", "6", "7", "8", "9"}
    assert set(bravais_summary["9"]).issuperset(
        {"bravais", "max_angular_difference", "unit_cell", "rmsd", "nspots"}
    )

    assert bravais_summary["9"]["unit_cell"] == pytest.approx(
        [57.78, 57.78, 150.0, 90.0, 90.0, 90.0], abs=1e-1
    )
    assert bravais_summary["9"]["bravais"] == "tP"
    assert bravais_summary["9"]["recommended"] is True
    assert bravais_summary["9"]["rmsd"] == pytest.approx(0.047, abs=1e-2)


def test_refine_bravais_settings_2(dials_regression, tmpdir):
    data_dir = os.path.join(dials_regression, "indexing_test_data", "multi_sweep")
    pickle_path = os.path.join(data_dir, "indexed.pickle")
    experiments_path = os.path.join(data_dir, "experiments.json")
    result = procrunner.run(
        ["dials.refine_bravais_settings", pickle_path, experiments_path],
        working_directory=tmpdir,
    )
    assert not result.returncode and not result.stderr
    for i in range(1, 10):
        assert tmpdir.join("bravais_setting_%i.expt" % i).check()

    experiments_list = load.experiment_list(
        tmpdir.join("bravais_setting_9.expt").strpath, check_format=False
    )
    assert len(experiments_list) == 4
    assert len(experiments_list.crystals()) == 1
    assert (
        experiments_list[0]
        .crystal.get_unit_cell()
        .is_similar_to(uctbx.unit_cell((7.31, 7.31, 6.82, 90.00, 90.00, 90.00)))
    )
    assert experiments_list[0].crystal.get_space_group().type().hall_symbol() == " I 4"
    assert tmpdir.join("bravais_summary.json").check()
    with tmpdir.join("bravais_summary.json").open("rb") as fh:
        bravais_summary = json.load(fh)
    for i in range(1, 23):
        assert str(i) in bravais_summary

    assert bravais_summary["9"]["unit_cell"] == pytest.approx(
        [7.31, 7.31, 6.82, 90.00, 90.00, 90.00], abs=1e-1
    )
    assert bravais_summary["9"]["bravais"] == "tI"
    assert bravais_summary["9"]["rmsd"] == pytest.approx(0.103, abs=1e-2)
    assert bravais_summary["9"]["recommended"] is True


def test_refine_bravais_settings_3(dials_regression, tmpdir):
    data_dir = os.path.join(dials_regression, "indexing_test_data", "trypsin")
    pickle_path = os.path.join(data_dir, "indexed.pickle")
    experiments_path = os.path.join(data_dir, "experiments.json")
    result = procrunner.run(
        [
            "dials.refine_bravais_settings",
            pickle_path,
            experiments_path,
            "crystal_id=1",
        ],
        working_directory=tmpdir,
    )
    assert not result.returncode and not result.stderr
    for i in range(1, 10):
        assert tmpdir.join("bravais_setting_%i.expt" % i).check()

    experiments_list = load.experiment_list(
        tmpdir.join("bravais_setting_5.expt").strpath, check_format=False
    )
    assert len(experiments_list) == 1
    assert (
        experiments_list[0]
        .crystal.get_unit_cell()
        .is_similar_to(uctbx.unit_cell((54.37, 58.29, 66.51, 90.00, 90.00, 90.00)))
    )
    assert (
        experiments_list[0].crystal.get_space_group().type().hall_symbol() == " P 2 2"
    )

    assert tmpdir.join("bravais_summary.json").check()
    with tmpdir.join("bravais_summary.json").open("rb") as fh:
        bravais_summary = json.load(fh)
    assert set(bravais_summary) == {"1", "2", "3", "4", "5", "6", "7", "8", "9"}

    assert bravais_summary["5"]["unit_cell"] == pytest.approx(
        [54.37, 58.29, 66.51, 90.00, 90.00, 90.00], abs=1e-1
    )
    assert bravais_summary["5"]["bravais"] == "oP"
    assert bravais_summary["5"]["rmsd"] == pytest.approx(0.1200, abs=1e-2)
    assert bravais_summary["5"]["recommended"] is True
    assert bravais_summary["9"]["recommended"] is False


def test_refine_bravais_settings_554(dials_regression, tmpdir):
    data_dir = os.path.join(dials_regression, "dials-554")
    pickle_path = os.path.join(data_dir, "indexed.pickle")
    experiments_path = os.path.join(data_dir, "experiments.json")
    result = procrunner.run(
        ["dials.refine_bravais_settings", pickle_path, experiments_path],
        working_directory=tmpdir,
    )
    assert not result.returncode and not result.stderr
    for i in range(1, 5):
        assert tmpdir.join("bravais_setting_%i.expt" % i).check()

    experiments_list = load.experiment_list(
        tmpdir.join("bravais_setting_5.expt").strpath, check_format=False
    )
    assert len(experiments_list) == 7
    assert len(experiments_list.crystals()) == 1
    crystal = experiments_list.crystals()[0]
    assert crystal.get_unit_cell().is_similar_to(
        uctbx.unit_cell((4.75863, 4.75863, 12.9885, 90, 90, 120))
    )
    assert crystal.get_space_group().type().hall_symbol() == " R 3"
    # assert all of the detectors are different
    for expt in experiments_list[1:]:
        assert expt.detector != experiments_list[0].detector
    for i in (0, 1, 6):
        assert experiments_list[i].detector[0].get_origin() == pytest.approx(
            (-41, 5.5, -135), abs=1
        )
    for i in (2, 3, 4, 5):
        assert experiments_list[i].detector[0].get_origin() == pytest.approx(
            (-41, 91, -99), abs=1
        )
    assert tmpdir.join("bravais_summary.json").check()
    with tmpdir.join("bravais_summary.json").open("rb") as fh:
        bravais_summary = json.load(fh)
    for i in range(1, 5):
        assert str(i) in bravais_summary

    assert bravais_summary["5"]["unit_cell"] == pytest.approx(
        [4.75863, 4.75863, 12.9885, 90, 90, 120], abs=1e-1
    )
    assert bravais_summary["5"]["bravais"] == "hR"
    assert bravais_summary["5"]["rmsd"] == pytest.approx(0.104, abs=1e-2)
    assert bravais_summary["5"]["recommended"] is True
