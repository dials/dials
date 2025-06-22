from __future__ import annotations

import shutil
import subprocess
from os import path

import pytest

from dxtbx.serialize import load

from dials.command_line.modify_experiments import phil_scope, update


def test_run(dials_data, tmp_path):
    orig_expt_json = dials_data("experiment_test_data") / "kappa_experiments.json"

    assert path.exists(orig_expt_json)

    orig_expt = load.experiment_list(orig_expt_json, check_format=False)
    orig_gonio = orig_expt.goniometers()[0]
    assert orig_gonio.get_angles() == pytest.approx([0, 180, 0])

    result = subprocess.run(
        [
            shutil.which("dials.modify_experiments"),
            orig_expt_json,
            "angles=10,20,30",
            "unit_cell=94,94,130,90,90,120",
        ],
        cwd=tmp_path,
        capture_output=True,
    )
    assert not result.returncode and not result.stderr

    new_expt_json = tmp_path / "modified.expt"
    assert new_expt_json.is_file()

    new_expt = load.experiment_list(new_expt_json, check_format=False)

    new_gonio = new_expt.goniometers()[0]
    assert new_gonio.get_angles() == pytest.approx([10, 20, 30])

    new_crystal = new_expt.crystals()[0]
    assert new_crystal.get_unit_cell().parameters() == pytest.approx(
        [94, 94, 130, 90, 90, 120]
    )


def test_select_experiments(dials_data, tmp_path):
    orig_expt_json = dials_data("indexing_test_data") / "multi_sweep-experiments.json"

    assert path.exists(orig_expt_json)

    result = subprocess.run(
        [
            shutil.which("dials.modify_experiments"),
            orig_expt_json,
            "select_experiments=0,1",
            "A_matrix=-0.076948,0.058256,0.104294,-0.010462,0.113451,-0.081650,-0.112936,-0.050201,-0.063496",
        ],
        cwd=tmp_path,
        capture_output=True,
    )
    assert not result.returncode and not result.stderr

    new_expt_json = tmp_path / "modified.expt"
    assert new_expt_json.is_file()

    new_expt = load.experiment_list(new_expt_json, check_format=False)

    # The original experiment list contains 4 experiments, all sharing a crystal model.
    # The selected experiments were copied so they do not share models. As a result,
    # there should now be 3 crystal models.
    assert len(new_expt.crystals()) == 3
    assert len(new_expt) == 4

    assert new_expt[0].crystal.get_A() == pytest.approx(
        [
            -0.076948,
            0.058256,
            0.104294,
            -0.010462,
            0.113451,
            -0.081650,
            -0.112936,
            -0.050201,
            -0.063496,
        ]
    )
    assert new_expt[1].crystal == new_expt[0].crystal

    assert new_expt[2].crystal != new_expt[0].crystal


def test_update(dials_data):
    orig_expt = dials_data("aluminium_standard", pathlib=True) / "imported.expt"
    assert orig_expt.is_file()

    orig_expt = load.experiment_list(orig_expt, check_format=False)
    orig_beam = orig_expt.beams()[0]
    assert orig_beam.get_wavelength() == pytest.approx(0.02508235604)

    working_params = phil_scope.fetch().extract()
    working_params.geometry.beam.wavelength = 0.05

    new_expt = update(orig_expt, working_params)

    new_beam = new_expt.beams()[0]
    assert new_beam.get_wavelength() == pytest.approx(0.05)
