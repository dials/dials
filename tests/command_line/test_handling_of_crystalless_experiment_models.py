from __future__ import annotations

import shutil
import subprocess

import pytest

from dxtbx.serialize import load

from dials.array_family import flex


@pytest.fixture
def indexed_data(dials_data, tmp_path):
    loc = dials_data("semisynthetic_multilattice", pathlib=True)
    result = subprocess.run(
        [
            shutil.which("dials.index"),
            loc / "ag_imported_1_50.expt",
            loc / "ag_strong_1_50.refl",
            loc / "bh_imported_1_50.expt",
            loc / "bh_strong_1_50.refl",
            "max_lattices=2",
            "joint_indexing=False",
            "n_macro_cycles=2",
            "output.retain_unindexed_experiments=True",
        ],
        cwd=tmp_path,
        capture_output=True,
    )
    assert not result.returncode and not result.stderr
    assert (tmp_path / "indexed.refl").is_file()
    assert (tmp_path / "indexed.expt").is_file()
    expts = load.experiment_list(tmp_path / "indexed.expt", check_format=False)
    assert len(expts) == 6
    assert len([c for c in expts.crystals() if c]) == 4
    refls = flex.reflection_table.from_file(tmp_path / "indexed.refl")
    refls.assert_experiment_identifiers_are_consistent(expts)
    assert set(refls["imageset_id"]) == {0, 1}
    return (tmp_path / "indexed.expt", tmp_path / "indexed.refl")


def test_refine(indexed_data, tmp_path):
    result = subprocess.run(
        [
            shutil.which("dials.refine"),
            indexed_data[0],
            indexed_data[1],
        ],
        cwd=tmp_path,
        capture_output=True,
    )
    assert not result.returncode and not result.stderr
    assert (tmp_path / "refined.refl").is_file()
    assert (tmp_path / "refined.expt").is_file()
    expts = load.experiment_list(tmp_path / "refined.expt", check_format=False)
    assert len(expts) == 6
    assert len([c for c in expts.crystals() if c]) == 4
    refls = flex.reflection_table.from_file(tmp_path / "refined.refl")
    refls.assert_experiment_identifiers_are_consistent(expts)
    assert set(refls["imageset_id"]) == {0, 1}

    result = subprocess.run(
        [
            shutil.which("dials.plot_scan_varying_model"),
            tmp_path / "refined.expt",
        ],
        cwd=tmp_path,
        capture_output=True,
    )
    assert not result.returncode and not result.stderr

    result = subprocess.run(
        [
            shutil.which("dials.filter_reflections"),
            tmp_path / "refined.expt",
            tmp_path / "refined.refl",
            "d_min=3.0",
        ],
        cwd=tmp_path,
        capture_output=True,
    )
    assert not result.returncode and not result.stderr
    refls = flex.reflection_table.from_file(tmp_path / "refined.refl")
    assert set(refls["id"]) == {0, 1, 2, 3, 4, 5}


@pytest.mark.parametrize(
    "program,requires_refl,params",
    [
        ("dials.stereographic_projection", False, ["hkl=1,0,0"]),
        ("dials.two_theta_offset", False, []),
        ("dials.unit_cell_histogram", False, []),
        ("dials.align_crystal", False, []),
        ("dials.cluster_unit_cell", False, []),
        ("dials.compare_orientation_matrices", False, []),
        ("dials.refine_bravais_settings", True, ["crystal_id=0"]),
        ("dials.indexed_as_integrated", True, []),
        ("dials.reindex", True, ["change_of_basis_op=a,b,c"]),
        ("dials.report", True, []),
        ("dials.show", True, []),
        ("dials.split_experiments", True, []),
        ("dials.spot_resolution_shells", True, []),
    ],
)
def test_programs(indexed_data, tmp_path, program, requires_refl, params):
    cmd = [shutil.which(program), indexed_data[0]]
    if requires_refl:
        cmd.append(indexed_data[1])
    if params:
        cmd.extend(params)
    result = subprocess.run(cmd, cwd=tmp_path, capture_output=True)
    assert not result.returncode and not result.stderr
