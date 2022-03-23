from __future__ import annotations

import os

import procrunner
import pytest

from dxtbx.serialize import load


def test_modify_geometry_run(dials_regression, tmp_path):
    orig_expt_json = os.path.join(
        dials_regression, "experiment_test_data/kappa_experiments.json"
    )
    assert orig_expt_json.is_file()
    orig_expt = load.experiment_list(orig_expt_json, check_format=False)

    result = procrunner.run(
        ["dials.modify_geometry", orig_expt_json, "angles=10,20,30"],
        working_directory=tmp_path,
    )
    assert not result.returncode and not result.stderr

    new_expt_json = tmp_path / "modified.expt"
    assert new_expt_json.is_file()

    new_expt = load.experiment_list(new_expt_json, check_format=False)

    orig_gonio = orig_expt.goniometers()[0]
    new_gonio = new_expt.goniometers()[0]
    assert orig_gonio.get_angles() == pytest.approx([0, 180, 0])
    assert new_gonio.get_angles() == pytest.approx([10, 20, 30])
