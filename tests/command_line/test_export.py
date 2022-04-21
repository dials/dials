from __future__ import annotations

import filecmp
import json

import procrunner
import pytest

import iotbx.cif
from cctbx import sgtbx, uctbx
from dxtbx.model.experiment_list import ExperimentListFactory
from dxtbx.serialize import load
from iotbx import mtz

from dials.array_family import flex
from dials.command_line.slice_sequence import slice_experiments, slice_reflections
from dials.util.multi_dataset_handling import assign_unique_identifiers


def run_export(export_format, dials_data, tmp_path):
    result = procrunner.run(
        [
            "dials.export",
            "format=" + export_format,
            dials_data("centroid_test_data", pathlib=True) / "experiments.json",
            dials_data("centroid_test_data", pathlib=True) / "integrated.pickle",
        ],
        working_directory=tmp_path,
    )
    assert not result.returncode and not result.stderr
    assert (tmp_path / f"integrated.{export_format}").is_file()


def test_nxs(dials_data, tmp_path):
    run_export("nxs", dials_data, tmp_path)


def test_mtz(dials_data, tmp_path):
    result = procrunner.run(
        [
            "dials.export",
            "format=mtz",
            "project_name=ham",
            "crystal_name=spam",
            dials_data("centroid_test_data", pathlib=True) / "experiments.json",
            dials_data("centroid_test_data", pathlib=True) / "integrated.pickle",
        ],
        working_directory=tmp_path,
    )
    assert not result.returncode and not result.stderr
    assert (tmp_path / "integrated.mtz").is_file()
    mtz_obj = mtz.object(str(tmp_path / "integrated.mtz"))
    assert mtz_obj.crystals()[1].name() == "spam"
    assert mtz_obj.crystals()[1].project_name() == "ham"


def test_mtz_recalculated_cell(dials_data, tmp_path):
    # First run dials.two_theta_refine to ensure that the crystals have
    # recalculated_unit_cell set
    scaled_expt = (
        dials_data("x4wide_processed", pathlib=True) / "AUTOMATIC_DEFAULT_scaled.expt"
    )
    scaled_refl = (
        dials_data("x4wide_processed", pathlib=True) / "AUTOMATIC_DEFAULT_scaled.refl"
    )
    result = procrunner.run(
        ["dials.two_theta_refine", scaled_expt, scaled_refl],
        working_directory=tmp_path,
    )
    assert (tmp_path / "refined_cell.expt").is_file()
    refined_expt = load.experiment_list(
        tmp_path / "refined_cell.expt", check_format=False
    )
    ttr_cell = refined_expt.crystals()[0].get_recalculated_unit_cell()

    d_min = 1.3
    result = procrunner.run(
        [
            "dials.export",
            "format=mtz",
            tmp_path / "refined_cell.expt",
            scaled_refl,
            f"d_min={d_min:f}",
        ],
        working_directory=tmp_path,
    )
    assert not result.returncode and not result.stderr
    assert (tmp_path / "scaled.mtz").is_file()
    # The resulting mtz should have the same unit cell set as the recalculated_unit_cell
    # from dials.two_theta_refine
    for ma in mtz.object(str(tmp_path / "scaled.mtz")).as_miller_arrays():
        assert ttr_cell.parameters() == pytest.approx(ma.unit_cell().parameters())
        assert ma.d_min() >= d_min


def test_mtz_best_unit_cell(dials_data, tmp_path):
    scaled_expt = (
        dials_data("x4wide_processed", pathlib=True) / "AUTOMATIC_DEFAULT_scaled.expt"
    )
    scaled_refl = (
        dials_data("x4wide_processed", pathlib=True) / "AUTOMATIC_DEFAULT_scaled.refl"
    )
    best_unit_cell = uctbx.unit_cell((42, 42, 39, 90, 90, 90))
    d_min = 1.5
    result = procrunner.run(
        [
            "dials.export",
            "format=mtz",
            scaled_expt,
            scaled_refl,
            f"d_min={d_min:f}",
            "best_unit_cell=%g,%g,%g,%g,%g,%g" % best_unit_cell.parameters(),
        ],
        working_directory=tmp_path,
    )
    assert not result.returncode and not result.stderr
    assert (tmp_path / "scaled.mtz").is_file()
    # The resulting mtz should have the best_unit_cell as input to dials.export
    for ma in mtz.object(str(tmp_path / "scaled.mtz")).as_miller_arrays():
        assert best_unit_cell.parameters() == pytest.approx(ma.unit_cell().parameters())
        assert ma.d_min() >= d_min


def test_multi_sequence_integrated_mtz(dials_data, tmp_path):
    """Test dials.export on multi-sequence integrated data."""
    # first combine two integrated files
    data = dials_data("multi_crystal_proteinase_k", pathlib=True)
    result = procrunner.run(
        [
            "dials.combine_experiments",
            data / "experiments_1.json",
            data / "reflections_1.pickle",
            data / "experiments_2.json",
            data / "reflections_2.pickle",
        ],
        working_directory=tmp_path,
    )

    assert not result.returncode and not result.stderr
    assert (tmp_path / "combined.refl").is_file()
    assert (tmp_path / "combined.expt").is_file()

    # now export
    result = procrunner.run(
        [
            "dials.export",
            "format=mtz",
            "mtz.hklout=integrated.mtz",
            tmp_path / "combined.refl",
            tmp_path / "combined.expt",
        ],
        working_directory=tmp_path,
    )
    assert not result.returncode and not result.stderr
    assert (tmp_path / "integrated.mtz").is_file()


def test_mtz_multi_wavelength(dials_data, tmp_path):
    """Test multi-wavelength mtz export"""
    # First make suitable input - multi datasets experiment list and reflection
    # table with different wavelengths
    mcp = dials_data("multi_crystal_proteinase_k", pathlib=True)
    exp_1 = load.experiment_list(mcp / "experiments_1.json", check_format=False)
    exp_2 = load.experiment_list(mcp / "experiments_2.json", check_format=False)
    refl_1 = flex.reflection_table.from_file(mcp / "reflections_1.pickle")
    refl_2 = flex.reflection_table.from_file(mcp / "reflections_2.pickle")

    exp_1[0].beam.set_wavelength(0.5)
    exp_2[0].beam.set_wavelength(1.0)

    exp_1.extend(exp_2)
    reflection_list = [refl_1, refl_2]
    exps, refls = assign_unique_identifiers(exp_1, reflection_list)
    joint_refl = flex.reflection_table()
    for r in refls:
        joint_refl.extend(r)
    exps.as_json(tmp_path / "tmp_exp.expt")
    joint_refl.as_file(tmp_path / "tmp_refl.refl")

    # Now run
    result = procrunner.run(
        [
            "dials.export",
            tmp_path / "tmp_exp.expt",
            tmp_path / "tmp_refl.refl",
            "format=mtz",
            "mtz.hklout=unmerged.mtz",
        ],
        environment_override={"DIALS_EXPORT_DO_NOT_CHECK_FORMAT": "True"},
        working_directory=tmp_path,
    )
    assert not result.returncode and not result.stderr
    assert (tmp_path / "unmerged.mtz").is_file()

    # Inspect output
    m = mtz.object(str(tmp_path / "unmerged.mtz")).crystals()
    n_batches = []
    wavelengths = []
    for crystal in m:
        for dataset in crystal.datasets():
            wavelengths.append(dataset.wavelength())
            n_batches.append(dataset.n_batches())
    assert n_batches == [0, 25, 25]  # base, dataset1, dataset2
    assert wavelengths == [0, 0.5, 1.0]  # base, dataset1, dataset2


def test_mtz_primitive_cell(dials_data, tmp_path):
    scaled_expt = dials_data("insulin_processed", pathlib=True) / "scaled.expt"
    scaled_refl = dials_data("insulin_processed", pathlib=True) / "scaled.refl"

    # First reindex to the primitive setting
    expts = load.experiment_list(scaled_expt, check_format=False)
    cs = expts[0].crystal.get_crystal_symmetry()
    cb_op = cs.change_of_basis_op_to_primitive_setting()
    procrunner.run(
        [
            "dials.reindex",
            scaled_expt,
            scaled_refl,
            f'change_of_basis_op="{cb_op}"',
        ],
        working_directory=tmp_path,
    )

    # Now export the reindexed experiments/reflections
    procrunner.run(
        ["dials.export", tmp_path / "reindexed.expt", tmp_path / "reindexed.refl"],
        working_directory=tmp_path,
    )

    mtz_obj = mtz.object(str(tmp_path / "scaled.mtz"))
    cs_primitive = cs.change_basis(cb_op)
    assert mtz_obj.space_group() == cs_primitive.space_group()
    refl = flex.reflection_table.from_file(scaled_refl)
    refl = refl.select(~refl.get_flags(refl.flags.bad_for_scaling, all=False))
    for ma in mtz_obj.as_miller_arrays():
        assert ma.crystal_symmetry().is_similar_symmetry(cs_primitive)
        assert ma.d_max_min() == pytest.approx(
            (flex.max(refl["d"]), flex.min(refl["d"]))
        )


@pytest.mark.parametrize("compress", [None, "gz", "bz2", "xz"])
@pytest.mark.parametrize("hklout", [None, "my.cif"])
def test_mmcif(compress, hklout, dials_data, tmp_path):
    # Call dials.export after integration

    command = [
        "dials.export",
        "format=mmcif",
        dials_data("centroid_test_data", pathlib=True) / "experiments.json",
        dials_data("centroid_test_data", pathlib=True) / "integrated.pickle",
    ]
    if hklout is not None:
        command.append(f"mmcif.hklout={hklout}")
    else:
        hklout = "integrated.cif"
    if compress is not None:
        command.append(f"mmcif.compress={compress}")
        hklin = hklout + "." + compress
    else:
        hklin = hklout
    result = procrunner.run(command, working_directory=tmp_path)
    assert not result.returncode and not result.stderr
    assert (tmp_path / hklin).is_file()


@pytest.mark.parametrize("pdb_version", ["v5", "v5_next"])
def test_mmcif_on_scaled_data(dials_data, tmp_path, pdb_version):
    """Call dials.export format=mmcif after scaling"""
    scaled_expt = (
        dials_data("x4wide_processed", pathlib=True) / "AUTOMATIC_DEFAULT_scaled.expt"
    )
    scaled_refl = (
        dials_data("x4wide_processed", pathlib=True) / "AUTOMATIC_DEFAULT_scaled.refl"
    )
    command = [
        "dials.export",
        "format=mmcif",
        scaled_expt,
        scaled_refl,
        "mmcif.hklout=scaled.mmcif",
        "compress=None",
        f"pdb_version={pdb_version}",
    ]
    result = procrunner.run(command, working_directory=tmp_path)
    assert not result.returncode and not result.stderr
    assert (tmp_path / "scaled.mmcif").is_file()

    model = iotbx.cif.reader(file_path=str(tmp_path / "scaled.mmcif")).model()
    if pdb_version == "v5":
        assert "_pdbx_diffrn_data_section.id" not in model["dials"].keys()
    elif pdb_version == "v5_next":
        assert "_pdbx_diffrn_data_section.id" in model["dials"].keys()


def test_mmcif_p1_narrow_wedge(dials_data, tmp_path):
    """Call dials.export format=mmcif after scaling"""
    data_dir = dials_data("x4wide_processed", pathlib=True)

    refl = flex.reflection_table.from_file(data_dir / "AUTOMATIC_DEFAULT_scaled.refl")
    refl = slice_reflections(refl, [(1, 3)])
    refl.as_file(tmp_path / "p1_narrow.refl")

    expts = load.experiment_list(
        data_dir / "AUTOMATIC_DEFAULT_scaled.expt", check_format=False
    )
    expts = slice_experiments(expts, [(1, 3)])
    expts[0].crystal.set_space_group(sgtbx.space_group())
    expts.as_file(tmp_path / "p1_narrow.expt")

    command = [
        "dials.export",
        "format=mmcif",
        tmp_path / "p1_narrow.expt",
        tmp_path / "p1_narrow.refl",
        "mmcif.hklout=scaled.mmcif",
        "compress=None",
    ]
    result = procrunner.run(command, working_directory=tmp_path)
    assert not result.returncode and not result.stderr
    assert (tmp_path / "scaled.mmcif").is_file()

    model = iotbx.cif.reader(file_path=str(tmp_path / "scaled.mmcif")).model()
    assert model["dials"]["_reflns.pdbx_redundancy"] == "1.0"
    assert model["dials"]["_reflns.pdbx_CC_half"] == "0.0"


def test_xds_ascii(dials_data, tmp_path):
    # Call dials.export
    result = procrunner.run(
        [
            "dials.export",
            "intensity=sum",
            "format=xds_ascii",
            dials_data("centroid_test_data", pathlib=True) / "experiments.json",
            dials_data("centroid_test_data", pathlib=True) / "integrated.pickle",
        ],
        working_directory=tmp_path,
    )
    assert not result.returncode and not result.stderr
    assert (tmp_path / "DIALS.HKL").is_file()

    psi_values = {
        (-9, 7, -10): 153.430361,
        (-5, 11, -26): 175.559441,
        (-4, 23, 24): 129.468070,
        (2, 10, 20): 147.947274,
    }

    with (tmp_path / "DIALS.HKL").open() as fh:
        for record in fh:
            if record.startswith("!"):
                continue
            tokens = record.split()
            hkl = tuple(map(int, tokens[:3]))
            if hkl not in psi_values:
                continue
            psi = float(tokens[-1])
            assert psi == pytest.approx(psi_values[hkl], abs=0.1)


def test_sadabs(dials_data, tmp_path):
    # Call dials.export
    result = procrunner.run(
        [
            "dials.export",
            "intensity=sum",
            "mtz.partiality_threshold=0.99",
            "format=sadabs",
            dials_data("centroid_test_data", pathlib=True) / "experiments.json",
            dials_data("centroid_test_data", pathlib=True) / "integrated.pickle",
        ],
        working_directory=tmp_path,
    )
    assert not result.returncode and not result.stderr
    assert (tmp_path / "integrated.sad").is_file()

    direction_cosines = {
        (-9, 7, -10): (0.51253, -0.72107, 0.84696, -0.68476, -0.14130, -0.10561),
        (-5, 11, -26): (0.51310, -0.62895, 0.84711, -0.59223, -0.13830, -0.50366),
        (-4, 23, 24): (0.51308, -0.60578, 0.84711, -0.31416, -0.13840, 0.73099),
        (2, 10, 20): (0.51239, -0.46605, 0.84693, -0.61521, -0.14204, 0.63586),
    }

    with (tmp_path / "integrated.sad").open() as fh:
        for record in fh:
            record = record.replace("-", " -")
            tokens = record.split()
            hkl = tuple(map(int, tokens[:3]))
            cosines = tuple(map(float, tokens[6:12]))
            if hkl not in direction_cosines:
                continue
            assert cosines == pytest.approx(direction_cosines[hkl], abs=0.001)


def test_json(dials_data, tmp_path):
    # Call dials.export
    result = procrunner.run(
        [
            "dials.export",
            "format=json",
            dials_data("centroid_test_data", pathlib=True)
            / "imported_experiments.json",
            dials_data("centroid_test_data", pathlib=True) / "strong.pickle",
        ],
        working_directory=tmp_path,
    )
    assert not result.returncode and not result.stderr
    assert (tmp_path / "rlp.json").is_file()

    d = json.load((tmp_path / "rlp.json").open("rb"))
    assert set(d) == {"imageset_id", "experiments", "rlp", "experiment_id"}
    assert d["rlp"][:3] == [0.123413, 0.576679, 0.186326], d["rlp"][:3]
    assert d["imageset_id"][0] == 0
    assert d["experiment_id"][0] == 0
    experiments = ExperimentListFactory.from_dict(d["experiments"])
    imgset = experiments.imagesets()
    assert len(imgset) == 1


def test_json_shortened(dials_data, tmp_path):
    # Call dials.export
    result = procrunner.run(
        [
            "dials.export",
            "format=json",
            dials_data("centroid_test_data", pathlib=True) / "experiments.json",
            dials_data("centroid_test_data", pathlib=True) / "integrated.pickle",
            "json.filename=integrated.json",
            "n_digits=4",
            "compact=False",
        ],
        working_directory=tmp_path,
    )
    assert not result.returncode and not result.stderr
    assert (tmp_path / "integrated.json").is_file()

    d = json.load((tmp_path / "integrated.json").open("rb"))
    assert "imageset_id" in d
    assert "rlp" in d
    assert "experiment_id" in d
    assert d["rlp"][:3] == [-0.5975, -0.6141, 0.4702], d["rlp"][:3]
    assert d["imageset_id"][0] == 0
    assert d["experiment_id"][0] == 0


def test_shelx(dials_data, tmp_path):
    # Call dials.export
    result = procrunner.run(
        [
            "dials.export",
            "intensity=scale",
            "format=shelx",
            dials_data("l_cysteine_4_sweeps_scaled", pathlib=True)
            / "scaled_20_25.expt",
            dials_data("l_cysteine_4_sweeps_scaled", pathlib=True)
            / "scaled_20_25.refl",
        ],
        working_directory=tmp_path,
    )
    assert not result.returncode and not result.stderr
    assert (tmp_path / "dials.hkl").is_file()

    intensities_sigmas = {
        (4, 2, 4): (173.14, 15.39),
        (3, -3, -3): (324.13, 25.92),
        (4, 0, 2): (876.02, 69.34),
        (3, -2, -1): (463.11, 36.76),
    }

    with (tmp_path / "dials.hkl").open() as fh:
        max_intensity = -9999.0
        for record in fh:
            tokens = record.split()
            hkl = tuple(map(int, tokens[:3]))
            i_sigi = tuple(map(float, tokens[3:5]))
            if hkl not in intensities_sigmas:
                if i_sigi[0] > max_intensity:
                    max_intensity = i_sigi[0]
                continue
            assert i_sigi == pytest.approx(intensities_sigmas[hkl], abs=0.001)
        assert max_intensity == pytest.approx(9999.00, abs=0.001)


def test_shelx_ins(dials_data, tmp_path):
    # Call dials.export
    result = procrunner.run(
        [
            "dials.export",
            "intensity=scale",
            "format=shelx",
            dials_data("l_cysteine_4_sweeps_scaled", pathlib=True)
            / "scaled_20_25.expt",
            dials_data("l_cysteine_4_sweeps_scaled", pathlib=True)
            / "scaled_20_25.refl",
        ],
        working_directory=tmp_path,
    )
    assert not result.returncode and not result.stderr
    assert (tmp_path / "dials.ins").is_file()

    cell_esds = {
        "CELL": (5.4815, 8.2158, 12.1457, 90.000, 90.000, 90.000),
        "ZERR": (0.0005, 0.0007, 0.0011, 0.003, 0.004, 0.004),
    }

    with (tmp_path / "dials.ins").open() as fh:
        for line in fh:
            tokens = line.split()
            instruction = tokens[0]
            if instruction in cell_esds:
                result = tuple(map(float, tokens[2:8]))
                assert result == pytest.approx(cell_esds[instruction], abs=0.001)


def test_shelx_ins_best_unit_cell(dials_data, tmp_path):
    # Call dials.export
    result = procrunner.run(
        [
            "dials.export",
            "intensity=scale",
            "format=shelx",
            "best_unit_cell=5,8,12,90,90,90",
            dials_data("l_cysteine_4_sweeps_scaled", pathlib=True)
            / "scaled_20_25.expt",
            dials_data("l_cysteine_4_sweeps_scaled", pathlib=True)
            / "scaled_20_25.refl",
        ],
        working_directory=tmp_path,
    )
    assert not result.returncode and not result.stderr
    assert (tmp_path / "dials.ins").is_file()

    cell_esds = {
        "CELL": (5.0, 8.0, 12.0, 90.0, 90.0, 90.0),
    }

    with (tmp_path / "dials.ins").open() as fh:
        for line in fh:
            tokens = line.split()
            instruction = tokens[0]
            assert instruction != "ZERR"
            if instruction in cell_esds:
                result = tuple(map(float, tokens[2:8]))
                assert result == pytest.approx(cell_esds[instruction], abs=0.001)


def test_export_sum_or_profile_only(dials_data, tmp_path):
    expt = dials_data("insulin_processed", pathlib=True) / "integrated.expt"
    refl = dials_data("insulin_processed", pathlib=True) / "integrated.refl"

    for remove in "prf", "sum":
        removed = tmp_path / f"removed_{remove}.refl"
        data = flex.reflection_table.from_file(refl)
        del data[f"intensity.{remove}.value"]
        del data[f"intensity.{remove}.variance"]
        data.as_file(removed)
        result = procrunner.run(
            ["dials.export", expt, removed, f"mtz.hklout=removed_{remove}.mtz"],
            working_directory=tmp_path,
        )
        assert not result.returncode and not result.stderr
        assert (tmp_path / f"removed_{remove}.mtz").is_file()


@pytest.mark.parametrize("intensity_choice", ["profile", "sum"])
def test_pets(dials_data, tmp_path, intensity_choice):
    expt = dials_data("quartz_processed", pathlib=True) / "integrated.expt"
    refl = dials_data("quartz_processed", pathlib=True) / "integrated.refl"
    # Call dials.export
    result = procrunner.run(
        [
            "dials.export",
            "intensity=scale",
            "format=pets",
            "id=0",
            "step=1",
            "n_merged=2",
            "intensity=" + intensity_choice,
            "filename_prefix=" + intensity_choice,
            expt,
            refl,
        ],
        working_directory=tmp_path,
    )
    assert not result.returncode and not result.stderr
    output = tmp_path / (intensity_choice + ".cif_pets")

    if intensity_choice == "profile":
        reference = dials_data("quartz_processed", pathlib=True) / "dials_prf.cif_pets"
    else:
        reference = dials_data("quartz_processed", pathlib=True) / "dials_dyn.cif_pets"
    assert filecmp.cmp(output, reference)
