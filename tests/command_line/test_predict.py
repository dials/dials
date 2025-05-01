from __future__ import annotations

import shutil
import subprocess
from pathlib import Path

from dials.array_family import flex


def plausible(table):
    # Check the reflection IDs
    assert "id" in table
    assert "miller_index" in table
    assert "s1" in table
    assert "xyzcal.px" in table
    assert "xyzcal.mm" in table
    for row in table.rows():
        assert row["id"] == 0
    return True


def test_static_prediction(dials_data, tmp_path):
    result = subprocess.run(
        [
            shutil.which("dials.predict"),
            str(
                dials_data("misc_regression", pathlib=True)
                / "prediction-static-crystal.expt"
            ),
        ],
        cwd=tmp_path,
        capture_output=True,
    )
    assert not result.returncode and not result.stderr

    table = flex.reflection_table.from_file(tmp_path / "predicted.refl")
    assert len(table) == 1996
    assert plausible(table)


def test_scan_varying_prediction(dials_data, tmp_path):
    result = subprocess.run(
        [
            shutil.which("dials.predict"),
            str(
                dials_data("misc_regression", pathlib=True)
                / "prediction-varying-crystal.expt"
            ),
        ],
        cwd=tmp_path,
        capture_output=True,
    )
    assert not result.returncode and not result.stderr

    table = flex.reflection_table.from_file(tmp_path / "predicted.refl")
    assert len(table) == 1934
    assert plausible(table)


def test_force_static_prediction(dials_data, tmp_path):
    result = subprocess.run(
        [
            shutil.which("dials.predict"),
            str(
                dials_data("misc_regression", pathlib=True)
                / "prediction-varying-crystal.expt"
            ),
            "force_static=True",
        ],
        cwd=tmp_path,
    )
    assert not result.returncode and not result.stderr

    table = flex.reflection_table.from_file(tmp_path / "predicted.refl")
    assert len(table) == 1996
    assert plausible(table)


def test_experiment_parameters(dials_data: Path, tmp_path):
    result = subprocess.run(
        [
            shutil.which("dials.predict"),
            str(
                dials_data("misc_regression", pathlib=True)
                / "prediction-varying-crystal.expt"
            ),
        ],
        cwd=tmp_path,
    )
    assert not result.returncode and not result.stderr

    table = flex.reflection_table.from_file(tmp_path / "predicted.refl")
    assert table.experiment_identifiers()
