from __future__ import annotations

import shutil
import subprocess

from dials.array_family import flex


def test_import_integrate_hkl(dials_data, tmp_path):
    data = dials_data("centroid_test_data", pathlib=True)
    result = subprocess.run(
        [
            shutil.which("dials.import_xds"),
            "input.method=reflections",
            data / "INTEGRATE.HKL",
            data / "experiments.json",
        ],
        cwd=tmp_path,
        capture_output=True,
    )
    assert not result.returncode and not result.stderr

    table = flex.reflection_table.from_file(tmp_path / "integrate_hkl.refl")
    assert "miller_index" in table
    assert "id" in table
    assert "panel" in table
    assert "xyzcal.px" in table
    assert "xyzobs.px.value" in table
    assert "intensity.cor.value" in table
    assert "intensity.cor.variance" in table
    assert len(table) == 174911


def test_import_spot_xds(dials_data, tmp_path):
    result = subprocess.run(
        [
            shutil.which("dials.import_xds"),
            "input.method=reflections",
            dials_data("centroid_test_data", pathlib=True) / "SPOT.XDS",
        ],
        cwd=tmp_path,
        capture_output=True,
    )
    assert not result.returncode and not result.stderr

    table = flex.reflection_table.from_file(tmp_path / "spot_xds.refl")
    assert "miller_index" in table
    assert "id" in table
    assert "panel" in table
    assert "xyzobs.px.value" in table
    assert "intensity.sum.value" in table
    assert len(table) == 742


def test_import_spot_xds_with_filtering(dials_data, tmp_path):
    result = subprocess.run(
        [
            shutil.which("dials.import_xds"),
            "input.method=reflections",
            dials_data("centroid_test_data", pathlib=True) / "SPOT.XDS",
            "remove_invalid=True",
        ],
        cwd=tmp_path,
        capture_output=True,
    )
    assert not result.returncode and not result.stderr

    table = flex.reflection_table.from_file(tmp_path / "spot_xds.refl")
    assert "miller_index" in table
    assert "id" in table
    assert "panel" in table
    assert "xyzobs.px.value" in table
    assert "intensity.sum.value" in table
    assert len(table) == 664


def test_from_xds_files(dials_data, tmp_path):
    # Import from the image files
    result = subprocess.run(
        [
            shutil.which("dials.import_xds"),
            "input.method=experiment",
            "output.filename=import_xds.expt",
            dials_data("centroid_test_data", pathlib=True),
        ],
        cwd=tmp_path,
        capture_output=True,
    )
    assert not result.returncode and not result.stderr

    assert tmp_path.joinpath("import_xds.expt").is_file()
