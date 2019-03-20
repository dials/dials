from __future__ import absolute_import, division, print_function

import json
import procrunner
import pytest
from dxtbx.serialize.load import _decode_dict

# Tests used to check for h5py
# May need to add this again if lack of this check causes issues.


def test_nxs(dials_data, tmpdir):
    # Call dials.export
    result = procrunner.run(
        [
            "dials.export",
            "format=nxs",
            dials_data("centroid_test_data").join("experiments.json").strpath,
            dials_data("centroid_test_data").join("integrated.pickle").strpath,
        ],
        working_directory=tmpdir.strpath,
    )
    assert result["exitcode"] == 0
    assert result["stderr"] == ""
    assert tmpdir.join("integrated.nxs").check(file=1)


def test_mtz(dials_data, tmpdir):
    # Call dials.export
    result = procrunner.run(
        [
            "dials.export",
            "format=mtz",
            dials_data("centroid_test_data").join("experiments.json").strpath,
            dials_data("centroid_test_data").join("integrated.pickle").strpath,
        ],
        working_directory=tmpdir.strpath,
    )
    assert result["exitcode"] == 0
    assert result["stderr"] == ""
    assert tmpdir.join("integrated.mtz").check(file=1)


def test_mmcif(dials_data, tmpdir):
    # Call dials.export after integration
    result = procrunner.run(
        [
            "dials.export",
            "format=mmcif",
            dials_data("centroid_test_data").join("experiments.json").strpath,
            dials_data("centroid_test_data").join("integrated.pickle").strpath,
        ],
        working_directory=tmpdir.strpath,
    )
    assert result["exitcode"] == 0
    assert result["stderr"] == ""
    assert tmpdir.join("integrated.cif").check(file=1)

    # TODO include similar test for exporting scaled data in mmcif format


def test_xds_ascii(dials_data, tmpdir):
    # Call dials.export
    result = procrunner.run(
        [
            "dials.export",
            "intensity=sum",
            "format=xds_ascii",
            dials_data("centroid_test_data").join("experiments.json").strpath,
            dials_data("centroid_test_data").join("integrated.pickle").strpath,
        ],
        working_directory=tmpdir.strpath,
    )
    assert result["exitcode"] == 0
    assert result["stderr"] == ""
    assert tmpdir.join("DIALS.HKL").check(file=1)

    psi_values = {
        (-9, 7, -10): 153.430361,
        (-5, 11, -26): 175.559441,
        (-4, 23, 24): 129.468070,
        (2, 10, 20): 147.947274,
    }

    with tmpdir.join("DIALS.HKL").open() as fh:
        for record in fh:
            if record.startswith("!"):
                continue
            tokens = record.split()
            hkl = tuple(map(int, tokens[:3]))
            if not hkl in psi_values:
                continue
            psi = float(tokens[-1])
            assert psi == pytest.approx(psi_values[hkl], abs=0.1)


def test_sadabs(dials_data, tmpdir):
    # Call dials.export
    result = procrunner.run(
        [
            "dials.export",
            "intensity=sum",
            "mtz.partiality_threshold=0.99",
            "format=sadabs",
            dials_data("centroid_test_data").join("experiments.json").strpath,
            dials_data("centroid_test_data").join("integrated.pickle").strpath,
        ],
        working_directory=tmpdir.strpath,
    )
    assert result["exitcode"] == 0
    assert result["stderr"] == ""
    assert tmpdir.join("integrated.sad").check(file=1)

    direction_cosines = {
        (-9, 7, -10): (0.51253, -0.72107, 0.84696, -0.68476, -0.14130, -0.10561),
        (-5, 11, -26): (0.51310, -0.62895, 0.84711, -0.59223, -0.13830, -0.50366),
        (-4, 23, 24): (0.51308, -0.60578, 0.84711, -0.31416, -0.13840, 0.73099),
        (2, 10, 20): (0.51239, -0.46605, 0.84693, -0.61521, -0.14204, 0.63586),
    }

    with tmpdir.join("integrated.sad").open() as fh:
        for record in fh:
            record = record.replace("-", " -")
            tokens = record.split()
            hkl = tuple(map(int, tokens[:3]))
            cosines = tuple(map(float, tokens[6:12]))
            if not hkl in direction_cosines:
                continue
            assert cosines == pytest.approx(direction_cosines[hkl], abs=0.001)


def test_json(dials_data, tmpdir):
    # Call dials.export
    result = procrunner.run(
        [
            "dials.export",
            "format=json",
            dials_data("centroid_test_data").join("datablock.json").strpath,
            dials_data("centroid_test_data").join("strong.pickle").strpath,
        ],
        working_directory=tmpdir.strpath,
    )
    assert result["exitcode"] == 0
    assert result["stderr"] == ""
    assert tmpdir.join("rlp.json").check(file=1)

    from dxtbx.model.experiment_list import ExperimentListFactory

    with tmpdir.join("rlp.json").open("rb") as f:
        d = json.load(f, object_hook=_decode_dict)
    assert d.keys() == ["imageset_id", "experiments", "rlp", "experiment_id"], d.keys()
    assert d["rlp"][:3] == [0.123454, 0.57687, 0.186465], d["rlp"][:3]
    assert d["imageset_id"][0] == 0
    assert d["experiment_id"][0] == 0
    experiments = ExperimentListFactory.from_dict(d["experiments"])
    imgset = experiments.imagesets()
    assert len(imgset) == 1


def test_json_shortened(dials_data, tmpdir):
    # Call dials.export
    result = procrunner.run(
        [
            "dials.export",
            "format=json",
            dials_data("centroid_test_data").join("experiments.json").strpath,
            dials_data("centroid_test_data").join("integrated.pickle").strpath,
            "json.filename=integrated.json",
            "n_digits=4",
            "compact=False",
        ],
        working_directory=tmpdir.strpath,
    )
    assert result["exitcode"] == 0
    assert result["stderr"] == ""
    assert tmpdir.join("integrated.json").check(file=1)

    with tmpdir.join("integrated.json").open("rb") as f:
        d = json.load(f)
    assert "imageset_id" in d.keys()
    assert "rlp" in d.keys()
    assert "experiment_id" in d.keys()
    assert d["rlp"][:3] == [-0.5975, -0.6141, 0.4702], d["rlp"][:3]
    assert d["imageset_id"][0] == 0
    assert d["experiment_id"][0] == 0
