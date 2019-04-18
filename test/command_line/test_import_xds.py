from __future__ import absolute_import, division, print_function

import six.moves.cPickle as pickle
import procrunner


def test_import_integrate_hkl(dials_data, tmpdir):
    from dials.array_family import flex  # noqa: F401 import dependency

    result = procrunner.run(
        [
            "dials.import_xds",
            "input.method=reflections",
            dials_data("centroid_test_data").join("INTEGRATE.HKL").strpath,
            dials_data("centroid_test_data").join("experiments.json").strpath,
        ],
        working_directory=tmpdir.strpath,
    )
    assert result["exitcode"] == 0
    assert result["stderr"] == ""

    table = flex.reflection_table.from_msgpack_file("integrate_hkl.mpack")

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
            dials_data("centroid_test_data").join("SPOT.XDS").strpath,
        ],
        working_directory=tmpdir.strpath,
    )
    assert result["exitcode"] == 0
    assert result["stderr"] == ""

    table = flex.reflection_table.from_msgpack_file("spot_xds.mpack")

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
            dials_data("centroid_test_data").join("SPOT.XDS").strpath,
            "remove_invalid=True",
        ],
        working_directory=tmpdir.strpath,
    )
    assert result["exitcode"] == 0
    assert result["stderr"] == ""

    table = flex.reflection_table.from_msgpack_file("spot_xds.mpack")

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
            "output.filename=import_experiments.json",
            dials_data("centroid_test_data").strpath,
        ],
        working_directory=tmpdir.strpath,
    )
    assert result["exitcode"] == 0
    assert result["stderr"] == ""

    assert tmpdir.join("import_experiments.json").check(file=1)
