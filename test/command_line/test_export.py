from __future__ import absolute_import, division, print_function

import json
import os

import procrunner
import pytest
from cctbx import uctbx
from dials.array_family import flex
from dials.util.multi_dataset_handling import assign_unique_identifiers
from dxtbx.model import ExperimentList
from dxtbx.serialize.load import _decode_dict
from dxtbx.serialize import load
from iotbx import mtz


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
    result = procrunner.run(
        [
            "dials.export",
            "format=mtz",
            "project_name=ham",
            "crystal_name=spam",
            dials_data("centroid_test_data").join("experiments.json"),
            dials_data("centroid_test_data").join("integrated.pickle"),
        ],
        working_directory=tmpdir,
    )
    assert not result.returncode and not result.stderr
    assert tmpdir.join("integrated.mtz").check(file=1)
    mtz_obj = mtz.object(tmpdir.join("integrated.mtz").strpath)
    assert mtz_obj.crystals()[1].name() == "spam"
    assert mtz_obj.crystals()[1].project_name() == "ham"


def test_mtz_recalculated_cell(dials_data, tmpdir):
    # First run dials.two_theta_refine to ensure that the crystals have
    # recalculated_unit_cell set
    scaled_expt = dials_data("x4wide_processed").join("AUTOMATIC_DEFAULT_scaled.expt")
    scaled_refl = dials_data("x4wide_processed").join("AUTOMATIC_DEFAULT_scaled.refl")
    result = procrunner.run(
        ["dials.two_theta_refine", scaled_expt, scaled_refl],
        working_directory=tmpdir,
    )
    assert tmpdir.join("refined_cell.expt").check(file=1)
    refined_expt = load.experiment_list(
        tmpdir.join("refined_cell.expt").strpath, check_format=False
    )
    ttr_cell = refined_expt.crystals()[0].get_recalculated_unit_cell()

    d_min = 1.3
    result = procrunner.run(
        [
            "dials.export",
            "format=mtz",
            tmpdir.join("refined_cell.expt"),
            scaled_refl,
            "d_min=%f" % d_min,
        ],
        working_directory=tmpdir,
    )
    assert not result.returncode and not result.stderr
    assert tmpdir.join("scaled.mtz").check(file=1)
    # The resulting mtz should have the same unit cell set as the recalculated_unit_cell
    # from dials.two_theta_refine
    for ma in mtz.object(tmpdir.join("scaled.mtz").strpath).as_miller_arrays():
        assert ttr_cell.parameters() == pytest.approx(ma.unit_cell().parameters())
        assert ma.d_min() >= d_min


def test_mtz_best_unit_cell(dials_data, tmpdir):
    scaled_expt = dials_data("x4wide_processed").join("AUTOMATIC_DEFAULT_scaled.expt")
    scaled_refl = dials_data("x4wide_processed").join("AUTOMATIC_DEFAULT_scaled.refl")
    best_unit_cell = uctbx.unit_cell((42, 42, 39, 90, 90, 90))
    d_min = 1.5
    result = procrunner.run(
        [
            "dials.export",
            "format=mtz",
            scaled_expt,
            scaled_refl,
            "d_min=%f" % d_min,
            "best_unit_cell=%g,%g,%g,%g,%g,%g" % best_unit_cell.parameters(),
        ],
        working_directory=tmpdir,
    )
    assert not result.returncode and not result.stderr
    assert tmpdir.join("scaled.mtz").check(file=1)
    # The resulting mtz should have the best_unit_cell as input to dials.export
    for ma in mtz.object(tmpdir.join("scaled.mtz").strpath).as_miller_arrays():
        assert best_unit_cell.parameters() == pytest.approx(ma.unit_cell().parameters())
        assert ma.d_min() >= d_min


def test_multi_sequence_integrated_mtz(dials_data, tmpdir):
    """Test dials.export on multi-sequence integrated data."""
    # first combine two integrated files
    result = procrunner.run(
        [
            "dials.combine_experiments",
            dials_data("multi_crystal_proteinase_k").join("experiments_1.json"),
            dials_data("multi_crystal_proteinase_k").join("reflections_1.pickle"),
            dials_data("multi_crystal_proteinase_k").join("experiments_2.json"),
            dials_data("multi_crystal_proteinase_k").join("reflections_2.pickle"),
        ],
        working_directory=tmpdir,
    )

    assert not result.returncode and not result.stderr
    assert tmpdir.join("combined.refl").check(file=1)
    assert tmpdir.join("combined.expt").check(file=1)

    # now export
    result = procrunner.run(
        [
            "dials.export",
            "format=mtz",
            "mtz.hklout=integrated.mtz",
            tmpdir.join("combined.refl").strpath,
            tmpdir.join("combined.expt").strpath,
        ],
        working_directory=tmpdir,
    )
    assert not result.returncode and not result.stderr
    assert tmpdir.join("integrated.mtz").check(file=1)


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
    refl_1 = flex.reflection_table.from_file(mcp.join("reflections_1.pickle").strpath)
    refl_2 = flex.reflection_table.from_file(mcp.join("reflections_2.pickle").strpath)

    exp_1[0].beam.set_wavelength(0.5)
    exp_2[0].beam.set_wavelength(1.0)

    exp_1.extend(exp_2)
    reflection_list = [refl_1, refl_2]
    exps, refls = assign_unique_identifiers(exp_1, reflection_list)
    joint_refl = flex.reflection_table()
    for r in refls:
        joint_refl.extend(r)
    exps.as_json("tmp_exp.expt")
    joint_refl.as_file("tmp_refl.refl")

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


def test_mtz_primitive_cell(dials_data, tmpdir):
    scaled_expt = dials_data("insulin_processed") / "scaled.expt"
    scaled_refl = dials_data("insulin_processed") / "scaled.refl"

    # First reindex to the primitive setting
    expts = ExperimentList.from_file(scaled_expt.strpath, check_format=False)
    cs = expts[0].crystal.get_crystal_symmetry()
    cb_op = cs.change_of_basis_op_to_primitive_setting()
    procrunner.run(
        [
            "dials.reindex",
            scaled_expt.strpath,
            scaled_refl.strpath,
            'change_of_basis_op="%s"' % cb_op,
        ],
        working_directory=tmpdir.strpath,
    )

    # Now export the reindexed experiments/reflections
    procrunner.run(
        ["dials.export", "reindexed.expt", "reindexed.refl"],
        working_directory=tmpdir.strpath,
    )

    mtz_obj = mtz.object(os.path.join(tmpdir.strpath, "scaled.mtz"))
    cs_primitive = cs.change_basis(cb_op)
    assert mtz_obj.space_group() == cs_primitive.space_group()
    refl = flex.reflection_table.from_file(scaled_refl.strpath)
    refl = refl.select(~refl.get_flags(refl.flags.bad_for_scaling, all=False))
    for ma in mtz_obj.as_miller_arrays():
        assert ma.crystal_symmetry().is_similar_symmetry(cs_primitive)
        assert ma.d_max_min() == pytest.approx(
            (flex.max(refl["d"]), flex.min(refl["d"]))
        )


@pytest.mark.parametrize("hklout", [None, "my.cif"])
def test_mmcif(hklout, dials_data, tmpdir):
    # Call dials.export after integration

    command = [
        "dials.export",
        "format=mmcif",
        dials_data("centroid_test_data").join("experiments.json").strpath,
        dials_data("centroid_test_data").join("integrated.pickle").strpath,
    ]
    if hklout is not None:
        command.append("mmcif.hklout=%s" % hklout)
    result = procrunner.run(command, working_directory=tmpdir.strpath)
    assert not result.returncode and not result.stderr
    if hklout is not None:
        assert tmpdir.join(hklout).check(file=1)
    else:
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
    assert set(d) == {"imageset_id", "experiments", "rlp", "experiment_id"}
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
