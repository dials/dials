from __future__ import absolute_import, division, print_function

import pytest
import os
from libtbx import easy_run
from libtbx import easy_pickle
from cctbx import sgtbx
from dxtbx.serialize import load


def test_reindex(dials_regression, run_in_tmpdir):
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
    command = " ".join(commands)
    print(command)

    result = easy_run.fully_buffered(command=command).raise_if_errors()
    old_reflections = easy_pickle.load(pickle_path)
    assert os.path.exists("reindexed.refl")
    new_reflections = easy_pickle.load("reindexed.refl")
    old_experiments = load.experiment_list(experiments_path, check_format=False)
    assert os.path.exists("reindexed.expt")
    new_experiments = load.experiment_list("reindexed.expt", check_format=False)
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
        "change_of_basis_op=%s" % str(cb_op),
        "output.experiments=P4.expt",
    ]
    command = " ".join(commands)
    print(command)
    result = easy_run.fully_buffered(command=command).raise_if_errors()
    # apply one of the symops from the space group
    cb_op = sgtbx.change_of_basis_op("-x,-y,z")
    commands = [
        "dials.reindex",
        "P4.expt",
        "change_of_basis_op=%s" % str(cb_op),
        "output.experiments=P4_reindexed.expt",
    ]
    command = " ".join(commands)
    print(command)
    result = easy_run.fully_buffered(command=command).raise_if_errors()
    new_experiments1 = load.experiment_list("P4_reindexed.expt", check_format=False)
    assert new_experiments1[0].crystal.get_A() == pytest.approx(
        old_experiments[0].crystal.change_basis(cb_op).get_A()
    )
    #
    cb_op = sgtbx.change_of_basis_op("-x,-y,z")
    commands = [
        "dials.reindex",
        "P4.expt",
        "change_of_basis_op=auto",
        "reference.experiments=P4_reindexed.expt",
        "output.experiments=P4_reindexed2.expt",
    ]
    command = " ".join(commands)
    print(command)
    result = easy_run.fully_buffered(command=command).raise_if_errors()
    new_experiments2 = load.experiment_list("P4_reindexed2.expt", check_format=False)
    assert new_experiments1[0].crystal.get_A() == pytest.approx(
        new_experiments2[0].crystal.get_A()
    )


def test_reindex_multi_sweep(dials_regression, run_in_tmpdir):
    data_dir = os.path.join(dials_regression, "indexing_test_data", "multi_sweep")
    pickle_path = os.path.join(data_dir, "indexed.pickle")
    experiments_path = os.path.join(data_dir, "experiments.json")
    commands = [
        "dials.reindex",
        pickle_path,
        experiments_path,
        "change_of_basis_op=x+y,x-z,y-z",
    ]
    command = " ".join(commands)
    print(command)

    result = easy_run.fully_buffered(command=command).raise_if_errors()
    old_reflections = easy_pickle.load(pickle_path)
    assert os.path.exists("reindexed.refl")
    new_reflections = easy_pickle.load("reindexed.refl")
    old_experiments = load.experiment_list(experiments_path, check_format=False)
    assert os.path.exists("reindexed.expt")
    new_experiments = load.experiment_list("reindexed.expt", check_format=False)
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
    tmpdir.chdir()

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
    command = " ".join(commands)
    print(command)

    _ = easy_run.fully_buffered(command=command).raise_if_errors()
    assert os.path.exists("P4.refl")
    assert os.path.exists("P4.expt")
    new_experiments = load.experiment_list("P4.expt", check_format=False)
    assert new_experiments[0].crystal.get_space_group().type().hall_symbol() == " P 4"

    # Now have something in P4, get another dataset in a different indexing scheme

    cb_op = sgtbx.change_of_basis_op("a,-b,-c")
    commands = [
        "dials.reindex",
        "P4.refl",
        "P4.expt",
        "change_of_basis_op=%s" % str(cb_op),
        "output.experiments=P4_reindexed.expt",
        "output.reflections=P4_reindexed.refl",
    ]
    command = " ".join(commands)
    print(command)
    _ = easy_run.fully_buffered(command=command).raise_if_errors()

    # now run reference reindexing
    commands = [
        "dials.reindex",
        "P4.refl",
        "P4.expt",
        "reference.experiments=P4_reindexed.expt",
        "reference.reflections=P4_reindexed.refl",
    ]
    command = " ".join(commands)
    print(command)
    _ = easy_run.fully_buffered(command=command).raise_if_errors()

    # expect reindexed_reflections to be same as P4_reindexed, not P4_reflections
    reindexed_reflections = easy_pickle.load("reindexed.refl")
    P4_reindexed = easy_pickle.load("P4_reindexed.refl")
    P4_reflections = easy_pickle.load("P4.refl")

    h1, k1, l1 = reindexed_reflections["miller_index"].as_vec3_double().parts()
    h2, k2, l2 = P4_reindexed["miller_index"].as_vec3_double().parts()
    h3, k3, l3 = P4_reflections["miller_index"].as_vec3_double().parts()

    # hkl1 and hkl2 should be same, as should have been reindexed by against the
    # reference, with the program determing a reindexing operator of a,-b,-c
    assert list(h1) == pytest.approx(list(h2))
    assert list(l1) == pytest.approx(list(l2))
    assert list(k1) == pytest.approx(list(k2))
    # h1 and h3 should be same, but not l and k, as these dataset should differ
    # by an a twinning operator of a,-b,-c
    assert list(h1) == pytest.approx(list(h3))
    assert list(l1) != pytest.approx(list(l3))
    assert list(k1) != pytest.approx(list(k3))
