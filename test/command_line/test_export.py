from __future__ import absolute_import, division, print_function

import json
import os

import procrunner
import pytest
from dials.array_family import flex
from dials.util.multi_dataset_handling import assign_unique_identifiers
from dxtbx.serialize.load import _decode_dict
from dxtbx.serialize import load
from iotbx import mtz

# Tests used to check for h5py
# May need to add this again if lack of this check causes issues.


def run_export(export_format, dials_data, tmpdir):
    result = procrunner.run(
        [
            "dials.export",
            "format=" + export_format,
            dials_data("centroid_test_data").join("experiments.json"),
            dials_data("centroid_test_data").join("integrated.pickle"),
        ],
        working_directory=tmpdir,
    )
    assert not result.returncode and not result.stderr
    assert tmpdir.join("integrated." + export_format).check(file=1)


def test_nxs(dials_data, tmpdir):
    run_export("nxs", dials_data, tmpdir)


def test_mtz(dials_data, tmpdir):
    run_export("mtz", dials_data, tmpdir)


def test_mtz_multi_wavelength(dials_data, run_in_tmpdir):
    """Test multi-wavelength mtz export"""
    # First make suitable input - multi datasets experiment list and reflection
    # table with different wavelengths
    mcp = dials_data("multi_crystal_proteinase_k")
    exp_1 = load.experiment_list(
        mcp.join("experiments_1.json").strpath, check_format=False
    )
    exp_2 = load.experiment_list(
        mcp.join("experiments_2.json").strpath, check_format=False
    )
    refl_1 = flex.reflection_table.from_pickle(mcp.join("reflections_1.pickle").strpath)
    refl_2 = flex.reflection_table.from_pickle(mcp.join("reflections_2.pickle").strpath)

    exp_1[0].beam.set_wavelength(0.5)
    exp_2[0].beam.set_wavelength(1.0)

    exp_1.extend(exp_2)
    reflection_list = [refl_1, refl_2]
    exps, refls = assign_unique_identifiers(exp_1, reflection_list)
    joint_refl = flex.reflection_table()
    for r in refls:
        joint_refl.extend(r)
    exps.as_json("tmp_exp.expt")
    joint_refl.as_pickle("tmp_refl.refl")

    # Now run
    result = procrunner.run(
        [
            "dials.export",
            "experiments=tmp_exp.expt",
            "reflections=tmp_refl.refl",
            "format=mtz",
            "mtz.hklout=unmerged.mtz",
        ],
        environment_override={"DIALS_EXPORT_DO_NOT_CHECK_FORMAT": "True"},
        working_directory=run_in_tmpdir.strpath,
    )
    assert not result.returncode and not result.stderr
    assert os.path.exists("unmerged.mtz")

    # Inspect output
    m = mtz.object("unmerged.mtz").crystals()
    n_batches = []
    wavelengths = []
    for crystal in m:
        for dataset in crystal.datasets():
            wavelengths.append(dataset.wavelength())
            n_batches.append(dataset.n_batches())
    assert n_batches == [0, 25, 25]  # base, dataset1, dataset2
    assert wavelengths == [0, 0.5, 1.0]  # base, dataset1, dataset2


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
    assert not result.returncode and not result.stderr
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
    assert not result.returncode and not result.stderr
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
            if hkl not in psi_values:
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
    assert not result.returncode and not result.stderr
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
            if hkl not in direction_cosines:
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
    assert not result.returncode and not result.stderr
    assert tmpdir.join("rlp.json").check(file=1)

    from dxtbx.model.experiment_list import ExperimentListFactory

    with tmpdir.join("rlp.json").open("rb") as f:
        d = json.load(f, object_hook=_decode_dict)
    assert list(d.keys()) == ["imageset_id", "experiments", "rlp", "experiment_id"]
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
    assert not result.returncode and not result.stderr
    assert tmpdir.join("integrated.json").check(file=1)

    with tmpdir.join("integrated.json").open("rb") as f:
        d = json.load(f)
    assert "imageset_id" in d
    assert "rlp" in d
    assert "experiment_id" in d
    assert d["rlp"][:3] == [-0.5975, -0.6141, 0.4702], d["rlp"][:3]
    assert d["imageset_id"][0] == 0
    assert d["experiment_id"][0] == 0
