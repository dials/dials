from __future__ import annotations

import json
import os
import shutil
import subprocess

import gemmi
import pytest

import iotbx.cif
from cctbx import sgtbx, uctbx
from dxtbx.model import Beam
from dxtbx.model.experiment_list import ExperimentListFactory
from dxtbx.serialize import load
from iotbx import mtz

from dials.array_family import flex
from dials.command_line.slice_sequence import slice_experiments, slice_reflections
from dials.util.multi_dataset_handling import assign_unique_identifiers


def run_export(export_format, dials_data, tmp_path):
    result = subprocess.run(
        [
            shutil.which("dials.export"),
            "format=" + export_format,
            dials_data("centroid_test_data") / "experiments.json",
            dials_data("centroid_test_data") / "integrated.pickle",
        ],
        cwd=tmp_path,
        capture_output=True,
    )
    assert not result.returncode and not result.stderr
    assert (tmp_path / f"integrated.{export_format}").is_file()


def test_nxs(dials_data, tmp_path):
    run_export("nxs", dials_data, tmp_path)


def test_mtz(dials_data, tmp_path):
    result = subprocess.run(
        [
            shutil.which("dials.export"),
            "format=mtz",
            "project_name=ham",
            "crystal_name=spam",
            dials_data("centroid_test_data") / "experiments.json",
            dials_data("centroid_test_data") / "integrated.pickle",
        ],
        cwd=tmp_path,
        capture_output=True,
    )
    assert not result.returncode and not result.stderr
    assert (tmp_path / "integrated.mtz").is_file()
    mtz_obj = mtz.object(str(tmp_path / "integrated.mtz"))
    assert mtz_obj.crystals()[1].name() == "spam"
    assert mtz_obj.crystals()[1].project_name() == "ham"


def test_mtz_recalculated_cell(dials_data, tmp_path):
    # First run dials.two_theta_refine to ensure that the crystals have
    # recalculated_unit_cell set
    scaled_expt = dials_data("x4wide_processed") / "AUTOMATIC_DEFAULT_scaled.expt"
    scaled_refl = dials_data("x4wide_processed") / "AUTOMATIC_DEFAULT_scaled.refl"
    result = subprocess.run(
        [shutil.which("dials.two_theta_refine"), scaled_expt, scaled_refl],
        cwd=tmp_path,
        capture_output=True,
    )
    assert (tmp_path / "refined_cell.expt").is_file()
    refined_expt = load.experiment_list(
        tmp_path / "refined_cell.expt", check_format=False
    )
    ttr_cell = refined_expt.crystals()[0].get_recalculated_unit_cell()

    d_min = 1.3
    result = subprocess.run(
        [
            shutil.which("dials.export"),
            "format=mtz",
            tmp_path / "refined_cell.expt",
            scaled_refl,
            f"d_min={d_min:f}",
        ],
        cwd=tmp_path,
        capture_output=True,
    )
    assert not result.returncode and not result.stderr
    assert (tmp_path / "scaled.mtz").is_file()
    # The resulting mtz should have the same unit cell set as the recalculated_unit_cell
    # from dials.two_theta_refine
    for ma in mtz.object(str(tmp_path / "scaled.mtz")).as_miller_arrays():
        assert ttr_cell.parameters() == pytest.approx(ma.unit_cell().parameters())
        assert ma.d_min() >= d_min


def test_mtz_best_unit_cell(dials_data, tmp_path):
    scaled_expt = dials_data("x4wide_processed") / "AUTOMATIC_DEFAULT_scaled.expt"
    scaled_refl = dials_data("x4wide_processed") / "AUTOMATIC_DEFAULT_scaled.refl"
    best_unit_cell = uctbx.unit_cell((42, 42, 39, 90, 90, 90))
    d_min = 1.5
    result = subprocess.run(
        [
            shutil.which("dials.export"),
            "format=mtz",
            scaled_expt,
            scaled_refl,
            f"d_min={d_min:f}",
            "best_unit_cell={:g},{:g},{:g},{:g},{:g},{:g}".format(
                *best_unit_cell.parameters()
            ),
        ],
        cwd=tmp_path,
        capture_output=True,
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
    data = dials_data("multi_crystal_proteinase_k")
    result = subprocess.run(
        [
            shutil.which("dials.combine_experiments"),
            data / "experiments_1.json",
            data / "reflections_1.pickle",
            data / "experiments_2.json",
            data / "reflections_2.pickle",
        ],
        cwd=tmp_path,
        capture_output=True,
    )

    assert not result.returncode and not result.stderr
    assert (tmp_path / "combined.refl").is_file()
    assert (tmp_path / "combined.expt").is_file()

    # now export
    result = subprocess.run(
        [
            shutil.which("dials.export"),
            "format=mtz",
            "mtz.hklout=integrated.mtz",
            tmp_path / "combined.refl",
            tmp_path / "combined.expt",
        ],
        cwd=tmp_path,
        capture_output=True,
    )
    assert not result.returncode and not result.stderr
    assert (tmp_path / "integrated.mtz").is_file()


def test_mtz_multi_wavelength(dials_data, tmp_path):
    """Test multi-wavelength mtz export"""
    # First make suitable input - multi datasets experiment list and reflection
    # table with different wavelengths
    mcp = dials_data("multi_crystal_proteinase_k")
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
    result = subprocess.run(
        [
            shutil.which("dials.export"),
            tmp_path / "tmp_exp.expt",
            tmp_path / "tmp_refl.refl",
            "format=mtz",
            "mtz.hklout=unmerged.mtz",
        ],
        env={"DIALS_EXPORT_DO_NOT_CHECK_FORMAT": "True", **os.environ},
        cwd=tmp_path,
        capture_output=True,
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
    scaled_expt = dials_data("insulin_processed") / "scaled.expt"
    scaled_refl = dials_data("insulin_processed") / "scaled.refl"

    # Export in I23
    subprocess.run(
        [
            shutil.which("dials.export"),
            scaled_expt,
            scaled_refl,
            "mtz.hklout=reference.mtz",
        ],
        cwd=tmp_path,
        capture_output=True,
    ).check_returncode()

    # Now reindex to the primitive setting
    expts = load.experiment_list(scaled_expt, check_format=False)
    cs = expts[0].crystal.get_crystal_symmetry()
    cb_op = cs.change_of_basis_op_to_primitive_setting()
    subprocess.run(
        [
            shutil.which("dials.reindex"),
            scaled_expt,
            scaled_refl,
            f'change_of_basis_op="{cb_op}"',
        ],
        cwd=tmp_path,
    ).check_returncode()

    # Export the reindexed experiments/reflections - internally this will
    # reindex back to I23
    subprocess.run(
        [
            shutil.which("dials.export"),
            tmp_path / "reindexed.expt",
            tmp_path / "reindexed.refl",
            "mtz.hklout=from_reindexed.mtz",
        ],
        cwd=tmp_path,
        capture_output=True,
    ).check_returncode()

    # Now read both files and check that relevant data are identical
    mtz1 = gemmi.read_mtz_file(str(tmp_path / "reference.mtz"))
    mtz2 = gemmi.read_mtz_file(str(tmp_path / "from_reindexed.mtz"))

    assert mtz1.spacegroup == mtz2.spacegroup
    assert mtz1.cell == mtz2.cell
    assert (mtz1.array == mtz2.array).all()
    for b1, b2 in zip(mtz1.batches, mtz2.batches):
        assert list(b1.ints) == list(b2.ints)
        assert list(b1.floats) == list(b2.floats)


@pytest.mark.parametrize("compress", [None, "gz", "bz2", "xz"])
@pytest.mark.parametrize("hklout", [None, "my.cif"])
def test_mmcif(compress, hklout, dials_data, tmp_path):
    # Call dials.export after integration

    command = [
        shutil.which("dials.export"),
        "format=mmcif",
        dials_data("centroid_test_data") / "experiments.json",
        dials_data("centroid_test_data") / "integrated.pickle",
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
    result = subprocess.run(command, cwd=tmp_path, capture_output=True)
    assert not result.returncode and not result.stderr
    assert (tmp_path / hklin).is_file()


@pytest.mark.skipif(
    shutil.which("gemmi") is None, reason="Could not find gemmi executable"
)
@pytest.mark.parametrize("pdb_version", ["v5", "v5_next"])
def test_mmcif_on_scaled_data(dials_data, tmp_path, pdb_version):
    """Call dials.export format=mmcif after scaling"""
    scaled_expt = dials_data("x4wide_processed") / "AUTOMATIC_DEFAULT_scaled.expt"
    scaled_refl = dials_data("x4wide_processed") / "AUTOMATIC_DEFAULT_scaled.refl"
    command = [
        shutil.which("dials.export"),
        "format=mmcif",
        scaled_expt,
        scaled_refl,
        "mmcif.hklout=scaled.mmcif",
        "compress=None",
        f"pdb_version={pdb_version}",
    ]
    result = subprocess.run(command, cwd=tmp_path, capture_output=True)
    assert not result.returncode and not result.stderr
    assert (tmp_path / "scaled.mmcif").is_file()

    model = iotbx.cif.reader(file_path=str(tmp_path / "scaled.mmcif")).model()
    if pdb_version == "v5":
        assert "_pdbx_diffrn_data_section.id" not in model["dials"].keys()
        # check that gemmi can understand the output
        cmd = [
            shutil.which("gemmi"),
            "cif2mtz",
            tmp_path / "scaled.mmcif",
            tmp_path / "test.mtz",
        ]
        result = subprocess.run(cmd, cwd=tmp_path, capture_output=True)
        assert not result.returncode and not result.stderr
        assert (tmp_path / "test.mtz").is_file()

    elif pdb_version == "v5_next":
        assert "_pdbx_diffrn_data_section.id" in model["dials"].keys()


def test_mmcif_p1_narrow_wedge(dials_data, tmp_path):
    """Call dials.export format=mmcif after scaling"""
    data_dir = dials_data("x4wide_processed")

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
        shutil.which("dials.export"),
        "format=mmcif",
        tmp_path / "p1_narrow.expt",
        tmp_path / "p1_narrow.refl",
        "mmcif.hklout=scaled.mmcif",
        "compress=None",
    ]
    result = subprocess.run(command, cwd=tmp_path, capture_output=True)
    assert not result.returncode and not result.stderr
    assert (tmp_path / "scaled.mmcif").is_file()

    model = iotbx.cif.reader(file_path=str(tmp_path / "scaled.mmcif")).model()
    assert model["dials"]["_reflns.pdbx_redundancy"] == "1.000"
    assert model["dials"]["_reflns.pdbx_CC_half"] == "0.0000"


def test_xds_ascii(dials_data, tmp_path):
    # Call dials.export
    result = subprocess.run(
        [
            shutil.which("dials.export"),
            "intensity=sum",
            "format=xds_ascii",
            dials_data("centroid_test_data") / "experiments.json",
            dials_data("centroid_test_data") / "integrated.pickle",
        ],
        cwd=tmp_path,
        capture_output=True,
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


def test_xds_ascii_scaled(dials_data, tmp_path):
    result = subprocess.run(
        [
            shutil.which("dials.export"),
            "format=xds_ascii",
            dials_data("l_cysteine_4_sweeps_scaled") / "scaled_30.expt",
            dials_data("l_cysteine_4_sweeps_scaled") / "scaled_30.refl",
        ],
        cwd=tmp_path,
        capture_output=True,
    )
    assert not result.returncode and not result.stderr
    assert (tmp_path / "DIALS.HKL").is_file()

    intensity_values = {
        (-7, 2, 11): 352.20 / 1.08375,
        (-8, 0, 2): 144.477 / 0.8753,
    }

    with (tmp_path / "DIALS.HKL").open() as fh:
        for record in fh:
            if record.startswith("!"):
                continue
            tokens = record.split()
            hkl = tuple(map(int, tokens[:3]))
            if hkl not in intensity_values:
                continue
            intensity = float(tokens[3])
            assert intensity == pytest.approx(intensity_values[hkl], abs=0.1)


def test_sadabs(dials_data, tmp_path):
    # Call dials.export
    result = subprocess.run(
        [
            shutil.which("dials.export"),
            "intensity=sum",
            "mtz.partiality_threshold=0.99",
            "format=sadabs",
            dials_data("centroid_test_data") / "experiments.json",
            dials_data("centroid_test_data") / "integrated.pickle",
        ],
        cwd=tmp_path,
        capture_output=True,
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
    result = subprocess.run(
        [
            shutil.which("dials.export"),
            "format=json",
            dials_data("centroid_test_data") / "imported_experiments.json",
            dials_data("centroid_test_data") / "strong.pickle",
        ],
        cwd=tmp_path,
        capture_output=True,
    )
    assert not result.returncode and not result.stderr
    assert (tmp_path / "rlp.json").is_file()

    with open(tmp_path / "rlp.json", mode="rb") as fh:
        d = json.load(fh)
    assert set(d) == {"imageset_id", "experiments", "rlp", "experiment_id"}
    assert d["rlp"][:3] == [0.123413, 0.576679, 0.186326], d["rlp"][:3]
    assert d["imageset_id"][0] == 0
    assert d["experiment_id"][0] == 0
    experiments = ExperimentListFactory.from_dict(d["experiments"])
    imgset = experiments.imagesets()
    assert len(imgset) == 1


def test_json_shortened(dials_data, tmp_path):
    # Call dials.export
    result = subprocess.run(
        [
            shutil.which("dials.export"),
            "format=json",
            dials_data("centroid_test_data") / "experiments.json",
            dials_data("centroid_test_data") / "integrated.pickle",
            "json.filename=integrated.json",
            "n_digits=4",
            "compact=False",
        ],
        cwd=tmp_path,
        capture_output=True,
    )
    assert not result.returncode and not result.stderr
    assert (tmp_path / "integrated.json").is_file()

    with open(tmp_path / "integrated.json", mode="rb") as fh:
        d = json.load(fh)
    assert "imageset_id" in d
    assert "rlp" in d
    assert "experiment_id" in d
    assert d["rlp"][:3] == [-0.5975, -0.6141, 0.4702], d["rlp"][:3]
    assert d["imageset_id"][0] == 0
    assert d["experiment_id"][0] == 0


def test_shelx(dials_data, tmp_path):
    # Call dials.export
    result = subprocess.run(
        [
            shutil.which("dials.export"),
            "intensity=scale",
            "format=shelx",
            dials_data("l_cysteine_4_sweeps_scaled") / "scaled_20_25.expt",
            dials_data("l_cysteine_4_sweeps_scaled") / "scaled_20_25.refl",
        ],
        cwd=tmp_path,
        capture_output=True,
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
    result = subprocess.run(
        [
            shutil.which("dials.export"),
            "intensity=scale",
            "format=shelx",
            dials_data("l_cysteine_4_sweeps_scaled") / "scaled_20_25.expt",
            dials_data("l_cysteine_4_sweeps_scaled") / "scaled_20_25.refl",
        ],
        cwd=tmp_path,
        capture_output=True,
    )
    assert not result.returncode and not result.stderr
    assert (tmp_path / "dials.ins").is_file()

    cell_esds = {
        "CELL": (5.48154, 8.21578, 12.14570, 90.0000, 90.0000, 90.0000),
        "ZERR": (0.00050, 0.00073, 0.00109, 0.0034, 0.0037, 0.0037),
    }

    with (tmp_path / "dials.ins").open() as fh:
        for line in fh:
            tokens = line.split()
            instruction = tokens[0]
            if instruction in cell_esds:
                result = tuple(map(float, tokens[2:8]))
                assert result == pytest.approx(cell_esds[instruction], abs=0.0001)


def test_shelx_ins_best_unit_cell(dials_data, tmp_path):
    # Call dials.export
    result = subprocess.run(
        [
            shutil.which("dials.export"),
            "intensity=scale",
            "format=shelx",
            "best_unit_cell=5,8,12,90,90,90",
            dials_data("l_cysteine_4_sweeps_scaled") / "scaled_20_25.expt",
            dials_data("l_cysteine_4_sweeps_scaled") / "scaled_20_25.refl",
        ],
        cwd=tmp_path,
        capture_output=True,
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


def test_shelx_ins_composition(dials_data, tmp_path):
    # Call dials.export
    result = subprocess.run(
        [
            shutil.which("dials.export"),
            "intensity=scale",
            "format=shelx",
            "composition=C3H7NO2S",
            dials_data("l_cysteine_4_sweeps_scaled") / "scaled_20_25.expt",
            dials_data("l_cysteine_4_sweeps_scaled") / "scaled_20_25.refl",
        ],
        cwd=tmp_path,
        capture_output=True,
    )
    assert not result.returncode and not result.stderr
    assert (tmp_path / "dials.ins").is_file()

    sfac_unit = {
        "SFAC": "C H N O S",
        "UNIT": "12 28 4 8 4",
    }

    with (tmp_path / "dials.ins").open() as fh:
        for line in fh:
            tokens = line.split()
            instruction = tokens[0]
            if instruction in sfac_unit:
                result = " ".join(tokens[1:6])
                assert result == sfac_unit[instruction]


def _export_cif(dials_data, tmp_path, *args, expt=None):
    data = dials_data("l_cysteine_4_sweeps_scaled")
    result = subprocess.run(
        [
            shutil.which("dials.export"),
            "intensity=scale",
            "format=cif",
            expt or data / "scaled_20_25.expt",
            data / "scaled_20_25.refl",
            *args,
        ],
        cwd=tmp_path,
        capture_output=True,
    )
    assert not result.returncode and not result.stderr
    assert (tmp_path / "dials.cif").is_file()
    return gemmi.cif.read_file(str(tmp_path / "dials.cif")).sole_block()


def test_cif(dials_data, tmp_path):
    block = _export_cif(dials_data, tmp_path)

    cell = {
        "_cell_length_a": "5.4815(5)",
        "_cell_length_b": "8.2158(7)",
        "_cell_length_c": "12.1457(11)",
        # angles are fixed by the symmetry of P 2 2 2, so carry no esds
        "_cell_angle_alpha": "90.0000",
        "_cell_angle_beta": "90.0000",
        "_cell_angle_gamma": "90.0000",
    }
    for name, value in cell.items():
        assert block.find_value(name) == value

    assert block.find_value("_space_group_IT_number") == "16"
    symops = list(block.find_loop("_space_group_symop_operation_xyz"))
    assert len(symops) == sgtbx.space_group_info("P 2 2 2").group().order_z()

    # Olex2 requires this exact set of column names, with intensity_u rather
    # than intensity_sigma
    table = block.find(
        "_diffrn_refln_",
        ["index_h", "index_k", "index_l", "intensity_net", "intensity_u"],
    )
    assert len(table) == int(block.find_value("_diffrn_reflns_number"))
    # a 0 0 0 index would make Olex2 discard the remaining reflections
    assert not any(row[0] == row[1] == row[2] == "0" for row in table)

    intensities_sigmas = {
        (4, 2, 4): (1695.9738, 150.7415),
        (3, -3, -3): (3174.9830, 253.9027),
        (4, 0, 2): (8580.8647, 679.1620),
        (3, -2, -1): (4536.3098, 360.0463),
    }
    max_intensity = 0.0
    for row in table:
        hkl = tuple(int(row[i]) for i in range(3))
        i_sigi = (float(row[3]), float(row[4]))
        max_intensity = max(max_intensity, i_sigi[0])
        if hkl in intensities_sigmas:
            assert i_sigi == pytest.approx(intensities_sigmas[hkl], abs=0.001)
    # unlike format=shelx, intensities are not rescaled to fit a fixed width
    assert max_intensity > 9999.0

    # no composition was given
    assert block.find_value("_chemical_formula_sum") is None
    # nor a scale group code
    assert not list(block.find_loop("_diffrn_refln_scale_group_code"))

    # the file is valid CIF as far as cctbx is concerned
    cs = iotbx.cif.builders.crystal_symmetry_builder(
        iotbx.cif.reader(file_path=str(tmp_path / "dials.cif")).model()["dials"]
    ).crystal_symmetry
    assert cs.space_group() == sgtbx.space_group_info("P 2 2 2").group()
    assert cs.unit_cell().parameters() == pytest.approx(
        (5.4815, 8.2158, 12.1457, 90, 90, 90), abs=0.001
    )


def test_cif_chemical_formula(dials_data, tmp_path):
    block = _export_cif(dials_data, tmp_path, "chemical_formula=C3H7NO2S")
    assert block.find_value("_chemical_formula_sum") == "'C3 H7 N O2 S'"
    assert float(block.find_value("_chemical_formula_weight")) == pytest.approx(
        121.16, abs=0.01
    )
    assert block.find_value("_cell_formula_units_Z") == "4"


def test_cif_scale_group_code(dials_data, tmp_path):
    block = _export_cif(dials_data, tmp_path, "cif.scale_group_code=True")
    codes = set(block.find_loop("_diffrn_refln_scale_group_code"))
    assert codes == {"1", "2"}


def test_cif_best_unit_cell(dials_data, tmp_path):
    block = _export_cif(dials_data, tmp_path, "best_unit_cell=5,8,12,90,90,90")
    # a supplied cell has no esds
    assert block.find_value("_cell_length_a") == "5.0000"
    assert block.find_value("_cell_length_b") == "8.0000"
    assert block.find_value("_cell_length_c") == "12.0000"


def test_cif_extra(dials_data, tmp_path):
    block = _export_cif(
        dials_data,
        tmp_path,
        "cif.extra=_diffrn_source_type=LaB6 gun",
        "cif.extra=_exptl_crystal_description=nanocrystal",
        "cif.extra=_exptl_crystal_colour=?",
        # this one overrides a value DIALS derived for itself
        "cif.extra=_diffrn_detector_type=ASI Timepix",
    )
    assert block.find_value("_diffrn_source_type") == "'LaB6 gun'"
    assert block.find_value("_exptl_crystal_description") == "nanocrystal"
    # the CIF special value for "unknown" is only meaningful unquoted
    assert block.find_value("_exptl_crystal_colour") == "?"
    assert block.find_value("_diffrn_detector_type") == "'ASI Timepix'"


def test_cif_extra_bad_syntax(dials_data, tmp_path):
    result = subprocess.run(
        [
            shutil.which("dials.export"),
            "intensity=scale",
            "format=cif",
            "cif.extra=no_leading_underscore=1",
            dials_data("l_cysteine_4_sweeps_scaled") / "scaled_20_25.expt",
            dials_data("l_cysteine_4_sweeps_scaled") / "scaled_20_25.refl",
        ],
        cwd=tmp_path,
        capture_output=True,
    )
    assert result.returncode
    assert b"Cannot interpret cif.extra" in result.stdout + result.stderr


def test_cif_electron_diffraction(dials_data, tmp_path):
    # There is no scaled electron diffraction data in dials_data, so relabel an
    # existing experiment list as if it had been collected with 200 kV electrons
    expts = ExperimentListFactory.from_json_file(
        str(dials_data("l_cysteine_4_sweeps_scaled") / "scaled_20_25.expt"),
        check_format=False,
    )
    for expt in expts:
        expt.beam.set_probe(Beam.get_probe_from_name("electron"))
        expt.beam.set_wavelength(0.0251)
    expts.as_file(tmp_path / "electron.expt")

    block = _export_cif(dials_data, tmp_path, expt=tmp_path / "electron.expt")

    assert block.find_value("_diffrn_radiation_probe") == "electron"
    assert block.find_value("_diffrn_radiation_wavelength") == "0.02510"
    # the accelerating voltage in kV, derived from the wavelength
    assert block.find_value("_diffrn_source_voltage") == "200"
    assert block.find_value("_diffrn_source") == "'electron gun'"
    assert (
        block.find_value("_diffrn_measurement_device")
        == "'transmission electron microscope'"
    )
    assert block.find_value("_diffrn_measurement_rotation_mode") == "rotation"
    assert (
        block.find_value("_diffrn_measurement_method")
        == "'continuous rotation 3D electron diffraction'"
    )


def test_export_sum_or_profile_only(dials_data, tmp_path):
    expt = dials_data("insulin_processed") / "integrated.expt"
    refl = dials_data("insulin_processed") / "integrated.refl"

    for remove in "prf", "sum":
        removed = tmp_path / f"removed_{remove}.refl"
        data = flex.reflection_table.from_file(refl)
        del data[f"intensity.{remove}.value"]
        del data[f"intensity.{remove}.variance"]
        data.as_file(removed)
        result = subprocess.run(
            [
                shutil.which("dials.export"),
                expt,
                removed,
                f"mtz.hklout=removed_{remove}.mtz",
            ],
            cwd=tmp_path,
            capture_output=True,
        )
        assert not result.returncode and not result.stderr
        assert (tmp_path / f"removed_{remove}.mtz").is_file()


@pytest.mark.parametrize("intensity_choice", ["profile", "sum"])
def test_pets(dials_data, tmp_path, intensity_choice):
    expt = dials_data("quartz_processed") / "integrated.expt"
    refl = dials_data("quartz_processed") / "integrated.refl"
    # Call dials.export
    result = subprocess.run(
        [
            shutil.which("dials.export"),
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
        cwd=tmp_path,
        capture_output=True,
    )
    assert not result.returncode and not result.stderr
    output = tmp_path / (intensity_choice + ".cif_pets")
    # On windows, \r\n endings are written; but our reference is \n
    output_data = output.read_bytes().replace(b"\r\n", b"\n")

    if intensity_choice == "profile":
        reference = dials_data("quartz_processed") / "dials_prf.cif_pets"
    else:
        reference = dials_data("quartz_processed") / "dials_dyn.cif_pets"
    assert output_data == reference.read_bytes()
