from __future__ import annotations

import shutil
import subprocess

from dials.array_family import flex


def test_import_integrate_hkl(dials_data, tmp_path):
    data = dials_data("centroid_test_data")
    result = subprocess.run(
        [shutil.which("dials.import_xds"), data / "INTEGRATE.HKL"],
        cwd=tmp_path,
        capture_output=True,
    )
    assert not result.returncode and not result.stderr
    assert (tmp_path / "integrate_hkl.refl").is_file()

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
    ## Make sure we can run this in a folder which just contains SPOT.XDS, as this
    ## is relied on by xia2 running XDS.
    shutil.copy(dials_data("centroid_test_data") / "SPOT.XDS", tmp_path)
    result = subprocess.run(
        [
            shutil.which("dials.import_xds"),
            tmp_path / "SPOT.XDS",
        ],
        cwd=tmp_path,
        capture_output=True,
    )
    assert not result.returncode and not result.stderr
    assert (tmp_path / "spot_xds.refl").is_file()
    assert not (
        tmp_path / "xds_models.expt"
    ).is_file()  # Don't expect this to be created in this case.

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
            dials_data("centroid_test_data") / "SPOT.XDS",
            "remove_invalid=True",
        ],
        cwd=tmp_path,
        capture_output=True,
    )
    assert not result.returncode and not result.stderr
    assert (tmp_path / "spot_xds.refl").is_file()

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
            dials_data("centroid_test_data"),
        ],
        cwd=tmp_path,
        capture_output=True,
    )
    assert not result.returncode and not result.stderr

    assert (tmp_path / "xds_models.expt").is_file()
    assert (tmp_path / "integrate_hkl.refl").is_file()


def test_from_xds_files_new(dials_data, tmp_path):
    # Import from more recent processed data.
    result = subprocess.run(
        [
            shutil.which("dials.import_xds"),
            dials_data("insulin_processed_xds_1"),
        ],
        cwd=tmp_path,
        capture_output=True,
    )
    assert not result.returncode and not result.stderr

    assert (tmp_path / "xds_models.expt").is_file()
    assert (tmp_path / "integrate_hkl.refl").is_file()

    # also use this as an opportunity to test mmcif export
    # of this data ends with the correct citation.
    result = subprocess.run(
        [
            shutil.which("dials.export"),
            tmp_path / "xds_models.expt",
            tmp_path / "integrate_hkl.refl",
            "format=mmcif",
        ],
        cwd=tmp_path,
        capture_output=True,
    )

    assert not result.returncode and not result.stderr
    assert (tmp_path / "integrated.cif").is_file()

    expected_phrase = "Data integration with the XDS package"
    with open(tmp_path / "integrated.cif", "r") as f:
        for line in f.readlines():
            if expected_phrase in line:
                return
        # reached end and wasn't found
        assert False, "Expected phrase not in integrated.cif"
