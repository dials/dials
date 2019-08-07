from __future__ import absolute_import, division, print_function

import pickle

import procrunner

from dials.array_family import flex


def test_import_integrate_hkl(dials_data, tmpdir):
    result = procrunner.run(
        [
            "dials.import_xds",
            "input.method=reflections",
            dials_data("centroid_test_data").join("INTEGRATE.HKL"),
            dials_data("centroid_test_data").join("experiments.json"),
        ],
        working_directory=tmpdir,
    )
    assert not result.returncode and not result.stderr

    with tmpdir.join("integrate_hkl.refl").open("rb") as fh:
        table = pickle.load(fh)

    assert "miller_index" in table
    assert "id" in table
    assert "panel" in table
    assert "xyzcal.px" in table
    assert "xyzobs.px.value" in table
    assert "intensity.cor.value" in table
    assert "intensity.cor.variance" in table
    assert len(table) == 174911


def test_import_spot_xds(dials_data, tmpdir):
    result = procrunner.run(
        [
            "dials.import_xds",
            "input.method=reflections",
            dials_data("centroid_test_data").join("SPOT.XDS"),
        ],
        working_directory=tmpdir,
    )
    assert not result.returncode and not result.stderr

    with tmpdir.join("spot_xds.refl").open("rb") as fh:
        table = pickle.load(fh)

    assert "miller_index" in table
    assert "id" in table
    assert "panel" in table
    assert "xyzobs.px.value" in table
    assert "intensity.sum.value" in table
    assert len(table) == 742


def test_import_spot_xds_with_filtering(dials_data, tmpdir):
    result = procrunner.run(
        [
            "dials.import_xds",
            "input.method=reflections",
            dials_data("centroid_test_data").join("SPOT.XDS"),
            "remove_invalid=True",
        ],
        working_directory=tmpdir,
    )
    assert not result.returncode and not result.stderr

    with tmpdir.join("spot_xds.refl").open("rb") as fh:
        table = pickle.load(fh)

    assert "miller_index" in table
    assert "id" in table
    assert "panel" in table
    assert "xyzobs.px.value" in table
    assert "intensity.sum.value" in table
    assert len(table) == 664


def test_from_xds_files(dials_data, tmpdir):
    # Import from the image files
    result = procrunner.run(
        [
            "dials.import_xds",
            "input.method=experiment",
            "output.filename=import_xds.expt",
            dials_data("centroid_test_data"),
        ],
        working_directory=tmpdir,
    )
    assert not result.returncode and not result.stderr

    assert tmpdir.join("import_xds.expt").check(file=1)
