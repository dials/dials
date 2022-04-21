from __future__ import annotations

import os

import procrunner
import pytest

import scitbx.matrix
from cctbx import sgtbx
from cctbx.sgtbx.lattice_symmetry import metric_subgroups
from dxtbx.model import Crystal, Experiment, ExperimentList
from dxtbx.serialize import load

import dials.command_line.reindex
from dials.array_family import flex
from dials.command_line.reindex import reindex_experiments


def test_reindex(dials_regression, tmpdir):
    data_dir = os.path.join(dials_regression, "indexing_test_data", "i04_weak_data")
    pickle_path = os.path.join(data_dir, "indexed.pickle")
    experiments_path = os.path.join(data_dir, "experiments.json")
    commands = [
        "dials.reindex",
        pickle_path,
        experiments_path,
        "change_of_basis_op=2a,b,c",
        "space_group=P1",
    ]
    result = procrunner.run(commands, working_directory=tmpdir)
    assert not result.returncode and not result.stderr

    old_reflections = flex.reflection_table.from_file(pickle_path)
    assert tmpdir.join("reindexed.refl").check()
    new_reflections = flex.reflection_table.from_file(
        tmpdir.join("reindexed.refl").strpath
    )
    old_experiments = load.experiment_list(experiments_path, check_format=False)
    assert tmpdir.join("reindexed.expt").check()
    new_experiments = load.experiment_list(
        tmpdir.join("reindexed.expt").strpath, check_format=False
    )
    h1, k1, l1 = old_reflections["miller_index"].as_vec3_double().parts()
    h2, k2, l2 = new_reflections["miller_index"].as_vec3_double().parts()
    assert 2 * h1 == pytest.approx(h2)
    assert k1 == pytest.approx(k2)
    assert l1 == pytest.approx(l2)
    old_uc_params = old_experiments[0].crystal.get_unit_cell().parameters()
    new_uc_params = new_experiments[0].crystal.get_unit_cell().parameters()
    assert new_uc_params[0] == pytest.approx(2 * old_uc_params[0])
    assert new_uc_params[1:] == pytest.approx(old_uc_params[1:])
    assert old_experiments[0].crystal.get_space_group().type().hall_symbol() == " P 1"
    assert new_experiments[0].crystal.get_space_group().type().hall_symbol() == " P 1"

    # set space group P4
    cb_op = sgtbx.change_of_basis_op("a,b,c")
    commands = [
        "dials.reindex",
        experiments_path,
        "space_group=P4",
        f"change_of_basis_op={cb_op}",
        "output.experiments=P4.expt",
    ]
    result = procrunner.run(commands, working_directory=tmpdir)
    assert not result.returncode and not result.stderr
    assert tmpdir.join("P4.expt").check()

    # apply one of the symops from the space group
    cb_op = sgtbx.change_of_basis_op("-x,-y,z")
    commands = [
        "dials.reindex",
        "P4.expt",
        f"change_of_basis_op={cb_op}",
        "output.experiments=P4_reindexed.expt",
    ]
    result = procrunner.run(commands, working_directory=tmpdir)
    assert not result.returncode and not result.stderr
    assert tmpdir.join("P4_reindexed.expt").check()

    new_experiments1 = load.experiment_list(
        tmpdir.join("P4_reindexed.expt").strpath, check_format=False
    )
    assert new_experiments1[0].crystal.get_A() == pytest.approx(
        old_experiments[0].crystal.change_basis(cb_op).get_A(), abs=1e-5
    )

    cb_op = sgtbx.change_of_basis_op("-x,-y,z")
    commands = [
        "dials.reindex",
        "P4.expt",
        "change_of_basis_op=auto",
        "reference.experiments=P4_reindexed.expt",
        "output.experiments=P4_reindexed2.expt",
    ]
    result = procrunner.run(commands, working_directory=tmpdir)
    assert not result.returncode and not result.stderr
    new_experiments2 = load.experiment_list(
        tmpdir.join("P4_reindexed2.expt").strpath, check_format=False
    )
    assert new_experiments1[0].crystal.get_A() == pytest.approx(
        new_experiments2[0].crystal.get_A()
    )

    # verify that the file can be read by dials.show
    commands = [
        "dials.show",
        tmpdir.join("reindexed.refl").strpath,
    ]
    result = procrunner.run(commands, working_directory=tmpdir)
    assert not result.returncode and not result.stderr


def test_reindex_multi_sequence(dials_regression, tmpdir):
    data_dir = os.path.join(dials_regression, "indexing_test_data", "multi_sweep")
    pickle_path = os.path.join(data_dir, "indexed.pickle")
    experiments_path = os.path.join(data_dir, "experiments.json")
    commands = [
        "dials.reindex",
        pickle_path,
        experiments_path,
        "change_of_basis_op=x+y,x-z,y-z",
    ]

    result = procrunner.run(commands, working_directory=tmpdir)
    assert not result.returncode and not result.stderr
    assert tmpdir.join("reindexed.refl").check()
    assert tmpdir.join("reindexed.expt").check()

    old_reflections = flex.reflection_table.from_file(pickle_path)
    new_reflections = flex.reflection_table.from_file(
        tmpdir.join("reindexed.refl").strpath
    )
    assert len(old_reflections) == len(new_reflections)
    new_experiments = load.experiment_list(
        tmpdir.join("reindexed.expt").strpath, check_format=False
    )
    new_cs = new_experiments[0].crystal.get_crystal_symmetry()
    assert new_cs.unit_cell().parameters() == pytest.approx(
        (
            6.189939294071243,
            6.189939294071243,
            6.189939294071242,
            113.16417286469935,
            107.65690626466579,
            107.65690626466579,
        )
    )
    assert (
        new_experiments[0].crystal.get_space_group().type().hall_symbol()
        == " I 4 (x+y,y+z,x+z)"
    )


def test_reindex_against_reference(dials_regression, tmpdir):
    """Test the reindexing against a reference dataset functionality."""
    data_dir = os.path.join(dials_regression, "indexing_test_data", "i04_weak_data")
    pickle_path = os.path.join(data_dir, "indexed.pickle")
    experiments_path = os.path.join(data_dir, "experiments.json")

    commands = [
        "dials.reindex",
        pickle_path,
        experiments_path,
        "change_of_basis_op=a,b,c",
        "space_group=P4",
        "output.reflections=P4.refl",
        "output.experiments=P4.expt",
    ]

    result = procrunner.run(commands, working_directory=tmpdir)
    assert not result.returncode and not result.stderr
    assert os.path.exists(tmpdir / "P4.refl")
    assert os.path.exists(tmpdir / "P4.expt")
    new_experiments = load.experiment_list(
        tmpdir.join("P4.expt").strpath, check_format=False
    )
    assert new_experiments[0].crystal.get_space_group().type().hall_symbol() == " P 4"

    # Now have something in P4, get another dataset in a different indexing scheme

    cb_op = sgtbx.change_of_basis_op("a,-b,-c")
    commands = [
        "dials.reindex",
        "P4.refl",
        "P4.expt",
        f"change_of_basis_op={cb_op}",
        "output.experiments=P4_reindexed.expt",
        "output.reflections=P4_reindexed.refl",
    ]
    result = procrunner.run(commands, working_directory=tmpdir)
    assert not result.returncode and not result.stderr

    # now run reference reindexing
    commands = [
        "dials.reindex",
        "P4.refl",
        "P4.expt",
        "reference.experiments=P4_reindexed.expt",
        "reference.reflections=P4_reindexed.refl",
    ]
    result = procrunner.run(commands, working_directory=tmpdir)
    assert not result.returncode and not result.stderr

    # expect reindexed_reflections to be same as P4_reindexed, not P4_reflections
    reindexed_reflections = flex.reflection_table.from_file(
        tmpdir.join("reindexed.refl").strpath
    )
    P4_reindexed = flex.reflection_table.from_file(
        tmpdir.join("P4_reindexed.refl").strpath
    )
    P4_reflections = flex.reflection_table.from_file(tmpdir.join("P4.refl").strpath)

    h1, k1, l1 = reindexed_reflections["miller_index"].as_vec3_double().parts()
    h2, k2, l2 = P4_reindexed["miller_index"].as_vec3_double().parts()
    h3, k3, l3 = P4_reflections["miller_index"].as_vec3_double().parts()

    # hkl1 and hkl2 should be same, as should have been reindexed by against the
    # reference, with the program determining a reindexing operator of a,-b,-c
    assert list(h1) == pytest.approx(list(h2))
    assert list(l1) == pytest.approx(list(l2))
    assert list(k1) == pytest.approx(list(k2))
    # h1 and h3 should be same, but not l and k, as these dataset should differ
    # by an a twinning operator of a,-b,-c
    assert list(h1) == pytest.approx(list(h3))
    assert list(l1) != pytest.approx(list(l3))
    assert list(k1) != pytest.approx(list(k3))


def test_reindex_experiments():
    # See also https://github.com/cctbx/cctbx_project/issues/424
    cs = sgtbx.space_group_info("I23").any_compatible_crystal_symmetry(volume=100000)
    B = scitbx.matrix.sqr(cs.unit_cell().fractionalization_matrix()).transpose()
    cryst = Crystal(B, cs.space_group())
    n_scan_points = 10
    A_at_scan_points = [(1, 0, 0, 0, 1, 0, 0, 0, 1)] * n_scan_points
    cryst.set_A_at_scan_points(A_at_scan_points)
    groups = metric_subgroups(cs, max_delta=5)
    for group in groups.result_groups:
        best_subsym = group["best_subsym"]
        cb_op = group["cb_op_inp_best"]
        expts = ExperimentList([Experiment(crystal=cryst)])
        reindexed_expts = reindex_experiments(
            experiments=expts, cb_op=cb_op, space_group=best_subsym.space_group()
        )
        assert (
            reindexed_expts[0]
            .crystal.get_crystal_symmetry()
            .is_similar_symmetry(best_subsym)
        )
        # Check that the scan-varying A matrices have been copied as well
        assert cryst.num_scan_points == n_scan_points


def test_reindex_cb_op_exit(dials_data):
    data_dir = dials_data("insulin_processed", pathlib=True)

    # Want a SystemExit, rather than an uncaught exception
    with pytest.raises(SystemExit) as e:
        dials.command_line.reindex.run(
            [
                "change_of_basis_op=a+2*b+c,-b+c,a",
                str(data_dir / "scaled.expt"),
                str(data_dir / "scaled.refl"),
            ]
        )
    # Make sure this is the expected error rather than e.g. an argument error
    assert "Unsuitable value" in str(e.value)


def test_reindex_reference_multi_crystal(dials_data, tmp_path):
    mcp = dials_data("multi_crystal_proteinase_k", pathlib=True)
    args = ["dials.cosym", "space_group=P4"]
    for i in [1, 2, 3, 4]:
        args.append(mcp / f"experiments_{i}.json")
        args.append(mcp / f"reflections_{i}.pickle")
    result = procrunner.run(args, working_directory=tmp_path)
    assert not result.returncode and not result.stderr
    assert tmp_path.joinpath("symmetrized.refl").is_file()
    assert tmp_path.joinpath("symmetrized.expt").is_file()

    args1 = [
        "dials.cosym",
        "space_group=P4",
        "output.experiments=sub_1.expt",
        "output.reflections=sub_1.refl",
    ]
    for i in [5, 7, 8, 10]:
        args1.append(mcp / f"experiments_{i}.json")
        args1.append(mcp / f"reflections_{i}.pickle")
    result1 = procrunner.run(args1, working_directory=tmp_path)
    assert not result1.returncode and not result1.stderr
    assert tmp_path.joinpath("sub_1.expt").is_file()
    assert tmp_path.joinpath("sub_1.refl").is_file()

    args2 = [
        "dials.reindex",
        "symmetrized.refl",
        "symmetrized.expt",
        "reference.experiments=sub_1.expt",
        "reference.reflections=sub_1.refl",
    ]
    result2 = procrunner.run(args2, working_directory=tmp_path)
    assert not result2.returncode and not result2.stderr
    assert tmp_path.joinpath("reindexed.expt").is_file()
    assert tmp_path.joinpath("reindexed.refl").is_file()
