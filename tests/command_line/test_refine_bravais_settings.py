from __future__ import annotations

import json
import os
import shutil
import subprocess
from pathlib import Path

import pytest

from cctbx import sgtbx, uctbx
from dxtbx.serialize import load

from dials.command_line import refine_bravais_settings


def test_refine_bravais_settings_i04_weak_data(dials_data, tmp_path):
    data_dir = dials_data("i04_weak_data")
    pickle_path = data_dir / "indexed.pickle"
    experiments_path = data_dir / "experiments.json"
    result = subprocess.run(
        [
            shutil.which("dials.refine_bravais_settings"),
            pickle_path,
            experiments_path,
            "reflections_per_degree=5",
            "minimum_sample_size=500",
            "beam.fix=all",
            "detector.fix=all",
            "prefix=tst_",
        ],
        cwd=tmp_path,
        capture_output=True,
    )
    assert not result.returncode and not result.stderr
    for i in range(1, 10):
        assert (tmp_path / f"tst_bravais_setting_{i}.expt").is_file()

    experiments_list = load.experiment_list(
        tmp_path / "tst_bravais_setting_9.expt", check_format=False
    )
    assert len(experiments_list) == 1
    assert (
        experiments_list[0]
        .crystal.get_unit_cell()
        .is_similar_to(uctbx.unit_cell((57.782, 57.782, 150.011, 90, 90, 90)))
    )
    assert experiments_list[0].crystal.get_space_group().type().hall_symbol() == " P 4"

    assert (tmp_path / "tst_bravais_summary.json").is_file()
    with (tmp_path / "tst_bravais_summary.json").open("rb") as fh:
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


def test_refine_bravais_settings_multi_sweep(dials_data, tmp_path):
    data_dir = dials_data("indexing_test_data")
    pickle_path = data_dir / "multi_sweep-indexed.pickle"
    experiments_path = data_dir / "multi_sweep-experiments.json"
    result = subprocess.run(
        [shutil.which("dials.refine_bravais_settings"), pickle_path, experiments_path],
        cwd=tmp_path,
        capture_output=True,
    )
    assert not result.returncode and not result.stderr
    for i in range(1, 10):
        assert (tmp_path / f"bravais_setting_{i}.expt").is_file()

    experiments_list = load.experiment_list(
        tmp_path / "bravais_setting_9.expt", check_format=False
    )
    assert len(experiments_list) == 4
    assert len(experiments_list.crystals()) == 1
    assert (
        experiments_list[0]
        .crystal.get_unit_cell()
        .is_similar_to(uctbx.unit_cell((7.31, 7.31, 6.82, 90.00, 90.00, 90.00)))
    )
    assert experiments_list[0].crystal.get_space_group().type().hall_symbol() == " I 4"
    assert (tmp_path / "bravais_summary.json").is_file()
    with (tmp_path / "bravais_summary.json").open("rb") as fh:
        bravais_summary = json.load(fh)
    for i in range(1, 23):
        assert str(i) in bravais_summary

    assert bravais_summary["9"]["unit_cell"] == pytest.approx(
        [7.31, 7.31, 6.82, 90.00, 90.00, 90.00], abs=1e-1
    )
    assert bravais_summary["9"]["bravais"] == "tI"
    assert bravais_summary["9"]["rmsd"] == pytest.approx(0.089, abs=1e-2)
    assert bravais_summary["9"]["recommended"] is True


def test_refine_bravais_settings_trypsin(dials_data: Path, tmp_path):
    data_dir = dials_data("indexing_test_data")
    pickle_path = os.path.join(data_dir, "trypsin-indexed.pickle")
    experiments_path = os.path.join(data_dir, "trypsin-experiments.json")
    result = subprocess.run(
        [
            shutil.which("dials.refine_bravais_settings"),
            pickle_path,
            experiments_path,
            "crystal_id=1",
        ],
        cwd=tmp_path,
        capture_output=True,
    )
    assert not result.returncode and not result.stderr
    for i in range(1, 10):
        assert (tmp_path / f"bravais_setting_{i}.expt").is_file()

    experiments_list = load.experiment_list(
        tmp_path / "bravais_setting_5.expt", check_format=False
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

    assert (tmp_path / "bravais_summary.json").is_file()
    with (tmp_path / "bravais_summary.json").open("rb") as fh:
        bravais_summary = json.load(fh)
    assert set(bravais_summary) == {"1", "2", "3", "4", "5", "6", "7", "8", "9"}

    assert bravais_summary["5"]["unit_cell"] == pytest.approx(
        [54.37, 58.29, 66.51, 90.00, 90.00, 90.00], abs=1e-1
    )
    assert bravais_summary["5"]["bravais"] == "oP"
    assert bravais_summary["5"]["rmsd"] == pytest.approx(0.103, abs=1e-2)
    assert bravais_summary["5"]["recommended"] is True
    assert bravais_summary["9"]["recommended"] is False


def test_refine_bravais_settings_554(dials_data, tmp_path):
    data_dir = dials_data("misc_regression", pathlib=True)
    reflections_path = str(data_dir / "dials-554_indexed.refl")
    experiments_path = str(data_dir / "dials-554_indexed.expt")
    result = subprocess.run(
        [
            shutil.which("dials.refine_bravais_settings"),
            reflections_path,
            experiments_path,
        ],
        cwd=tmp_path,
        capture_output=True,
    )
    assert not result.returncode and not result.stderr
    for i in range(1, 5):
        assert (tmp_path / f"bravais_setting_{i}.expt").is_file()

    experiments_list = load.experiment_list(
        tmp_path / "bravais_setting_5.expt", check_format=False
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
    assert (tmp_path / "bravais_summary.json").is_file()
    with (tmp_path / "bravais_summary.json").open("rb") as fh:
        bravais_summary = json.load(fh)
    for i in range(1, 5):
        assert str(i) in bravais_summary

    assert bravais_summary["5"]["unit_cell"] == pytest.approx(
        [4.75863, 4.75863, 12.9885, 90, 90, 120], abs=1e-1
    )
    assert bravais_summary["5"]["bravais"] == "hR"
    assert bravais_summary["5"]["rmsd"] == pytest.approx(0.104, abs=1e-2)
    assert bravais_summary["5"]["recommended"] is True


@pytest.mark.parametrize(
    "best_monoclinic_beta,expected_space_group,expected_unit_cell",
    [
        (True, "I 1 2 1", (44.47, 52.85, 111.46, 90.00, 99.91, 90.00)),
        (False, "C 1 2 1", (112.67, 52.85, 44.47, 90.00, 102.97, 90.00)),
    ],
)
def test_setting_c2_vs_i2(
    best_monoclinic_beta,
    expected_space_group,
    expected_unit_cell,
    dials_data,
    tmp_path,
):
    data_dir = dials_data("mpro_x0305_processed", pathlib=True)
    refl_path = data_dir / "indexed.refl"
    experiments_path = data_dir / "indexed.expt"
    result = subprocess.run(
        [
            shutil.which("dials.refine_bravais_settings"),
            os.fspath(experiments_path),
            os.fspath(refl_path),
            f"best_monoclinic_beta={best_monoclinic_beta}",
        ],
        cwd=tmp_path,
        capture_output=True,
    )
    assert not result.returncode and not result.stderr
    expts_orig = load.experiment_list(experiments_path, check_format=False)
    expts = load.experiment_list(
        tmp_path / "bravais_setting_2.expt", check_format=False
    )
    expected_bravais_lattice = str(
        sgtbx.bravais_types.bravais_lattice(symbol=expected_space_group)
    )
    expts[0].crystal.get_space_group().type().lookup_symbol() == expected_space_group
    assert expts[0].crystal.get_unit_cell().parameters() == pytest.approx(
        expected_unit_cell, abs=1e-2
    )
    with (tmp_path / "bravais_summary.json").open("rb") as fh:
        bravais_summary = json.load(fh)
    assert bravais_summary["2"]["bravais"] == expected_bravais_lattice
    # Verify that the cb_op converts from the input setting to the refined setting
    cb_op = sgtbx.change_of_basis_op(str(bravais_summary["2"]["cb_op"]))
    assert (
        expts_orig[0]
        .crystal.change_basis(cb_op)
        .get_unit_cell()
        .is_similar_to(
            expts[0].crystal.get_unit_cell(),
            relative_length_tolerance=0.1,
            absolute_angle_tolerance=1,
        )
    )
    stdout = result.stdout.decode()
    assert bravais_summary["2"]["cb_op"] in stdout
    assert f"| {expected_bravais_lattice}        |" in stdout
    expected_short_name = refine_bravais_settings.short_space_group_name(
        sgtbx.space_group_info(expected_space_group).group()
    )
    assert f"{expected_bravais_lattice}: {expected_short_name}" in stdout


def test_refine_bravais_settings_non_primitive_input(dials_data, tmp_path):
    data_dir = dials_data("insulin_processed", pathlib=True)
    refl_path = data_dir / "indexed.refl"
    expt_path = data_dir / "indexed.expt"
    result = subprocess.run(
        [
            shutil.which("dials.refine_bravais_settings"),
            expt_path,
            refl_path,
        ],
        cwd=tmp_path,
        capture_output=True,
    )
    assert not result.returncode and not result.stderr
    for i in range(1, 22):
        assert (tmp_path / f"bravais_setting_{i}.expt").is_file()

    assert (tmp_path / "bravais_summary.json").is_file()
    with (tmp_path / "bravais_summary.json").open("rb") as fh:
        bravais_summary = json.load(fh)
    assert bravais_summary["1"]["bravais"] == "aP"
    assert bravais_summary["1"]["cb_op"] == "y+z,x+z,x+y"
    assert bravais_summary["22"]["bravais"] == "cI"
    assert bravais_summary["22"]["cb_op"] == "a,b,c"

    for record in result.stdout.decode().split("\n"):
        if "aP " in record:
            assert "y+z,x+z,x+y" in record
        elif "cI " in record:
            assert "a,b,c" in record

    # Verify that the reported cb_ops correctly map the input setting
    # to the unit cell for each Bravais setting
    input_expts = load.experiment_list(expt_path, check_format=False)
    input_cs = input_expts[0].crystal.get_crystal_symmetry()
    for i in range(22):
        uc_input_to_ref = input_cs.unit_cell().change_basis(
            sgtbx.change_of_basis_op(bravais_summary[f"{i + 1}"]["cb_op"])
        )
        uc_ref = uctbx.unit_cell(bravais_summary[f"{i + 1}"]["unit_cell"])
        assert uc_input_to_ref.is_similar_to(uc_ref)
