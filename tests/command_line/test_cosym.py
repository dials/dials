from __future__ import annotations

import os
import pathlib

import pytest

from cctbx import sgtbx, uctbx
from dxtbx.serialize import load

import dials.command_line.cosym as dials_cosym
from dials.algorithms.symmetry.cosym._generate_test_data import (
    generate_experiments_reflections,
)
from dials.array_family import flex
from dials.util import Sorry


@pytest.mark.parametrize(
    "space_group,engine", [(None, "scitbx"), ("P 1", "scipy"), ("P 4", "scipy")]
)
def test_cosym(dials_data, run_in_tmp_path, space_group, engine):
    mcp = dials_data("multi_crystal_proteinase_k", pathlib=True)
    args = ["space_group=" + str(space_group), "seed=0", f"engine={engine}"]
    for i in [1, 2, 3, 4, 5, 7, 8, 10]:
        args.append(os.fspath(mcp / f"experiments_{i}.json"))
        args.append(os.fspath(mcp / f"reflections_{i}.pickle"))
    dials_cosym.run(args=args)
    assert pathlib.Path("symmetrized.refl").is_file()
    assert pathlib.Path("symmetrized.expt").is_file()
    experiments = load.experiment_list("symmetrized.expt", check_format=False)
    if space_group is None:
        assert (
            experiments[0].crystal.get_space_group().type().lookup_symbol() == "P 4 2 2"
        )
    else:
        assert (
            experiments[0].crystal.get_space_group().type().lookup_symbol()
            == space_group
        )
    joint_reflections = flex.reflection_table.from_file("symmetrized.refl")
    # check that there are 8 unique id and imageset_ids, and that these
    # correctly correspond to each experiment
    assert len(set(joint_reflections["id"])) == 8
    assert len(set(joint_reflections["imageset_id"])) == 8
    for id_ in range(8):
        sel = joint_reflections["id"] == id_
        assert set(joint_reflections["imageset_id"].select(sel)) == {id_}


def test_cosym_partial_dataset(dials_data, run_in_tmp_path):
    """Test how cosym handles partial/bad datasets."""
    mcp = dials_data("multi_crystal_proteinase_k", pathlib=True)
    args = []
    for i in [1, 2]:
        args.append(os.fspath(mcp / f"experiments_{i}.json"))
        args.append(os.fspath(mcp / f"reflections_{i}.pickle"))
    # Make one dataset that will be removed in prefiltering
    r = flex.reflection_table.from_file(mcp / "reflections_8.pickle")
    r["partiality"] = flex.double(r.size(), 0.1)
    r.as_file("renamed.refl")
    args.append("renamed.refl")
    args.append(os.fspath(mcp / "experiments_8.json"))
    # Add another good dataset at the end of the input list
    args.append(os.fspath(mcp / "experiments_10.json"))
    args.append(os.fspath(mcp / "reflections_10.pickle"))
    args.append("threshold=None")

    dials_cosym.run(args=args)
    assert pathlib.Path("symmetrized.refl").is_file()
    assert pathlib.Path("symmetrized.expt").is_file()
    experiments = load.experiment_list("symmetrized.expt", check_format=False)
    assert len(experiments) == 3


def test_cosym_resolution_filter_excluding_datasets(dials_data, run_in_tmp_path):
    mcp = dials_data("multi_crystal_proteinase_k", pathlib=True)
    args = ["space_group=P4", "seed=0", "d_min=20.0"]
    for i in [1, 2, 3, 4, 5, 7, 8, 10]:
        args.append(os.fspath(mcp / f"experiments_{i}.json"))
        args.append(os.fspath(mcp / f"reflections_{i}.pickle"))
    dials_cosym.run(args=args)
    assert pathlib.Path("symmetrized.refl").is_file()
    assert pathlib.Path("symmetrized.expt").is_file()

    joint_reflections = flex.reflection_table.from_file("symmetrized.refl")
    # check that there are 8 unique id and imageset_ids, and that these
    # correctly correspond to each experiment
    assert set(joint_reflections["id"]) == {0, 1, 2, 3, 4, 5, 6}
    assert len(set(joint_reflections["imageset_id"])) == 7
    for id_ in range(7):
        sel = joint_reflections["id"] == id_
        assert set(joint_reflections["imageset_id"].select(sel)) == {id_}


def test_cosym_partial_dataset_raises_sorry(dials_data, run_in_tmp_path, capsys):
    """Test how cosym handles partial/bad datasets."""
    mcp = dials_data("multi_crystal_proteinase_k", pathlib=True)
    args = ["renamed.refl", os.fspath(mcp / "experiments_8.json")]
    r2 = flex.reflection_table.from_file(mcp / "reflections_10.pickle")
    r2["partiality"] = flex.double(r2.size(), 0.1)
    r2.as_file("renamed2.refl")
    args.append("renamed2.refl")
    args.append(os.fspath(mcp / "experiments_10.json"))

    with pytest.raises(Sorry):
        dials_cosym.run(args=args)


@pytest.mark.parametrize(
    (
        "space_group",
        "unit_cell",
        "dimensions",
        "sample_size",
        "use_known_space_group",
        "use_known_lattice_group",
    ),
    [
        ("P2", None, None, 10, False, False),
        ("P3", None, None, 20, False, False),
        ("I23", None, 2, 10, False, False),
        ("P422", (79, 79, 37, 90, 90, 90), None, 10, True, False),
        ("P321", (59.39, 59.39, 28.35, 90, 90, 120), None, 5, False, False),
    ],
)
def test_synthetic(
    space_group,
    unit_cell,
    dimensions,
    sample_size,
    use_known_space_group,
    use_known_lattice_group,
    run_in_tmp_path,
):
    space_group = sgtbx.space_group_info(space_group).group()
    if unit_cell is not None:
        unit_cell = uctbx.unit_cell(unit_cell)
    experiments, reflections, _ = generate_experiments_reflections(
        space_group=space_group,
        unit_cell=unit_cell,
        unit_cell_volume=10000,
        sample_size=sample_size,
        map_to_p1=True,
        d_min=1.5,
    )

    experiments.as_json("tmp.expt")
    expt_file = "tmp.expt"
    joint_table = flex.reflection_table()
    for r in reflections:
        joint_table.extend(r)
    joint_table.as_file("tmp.refl")
    refl_file = "tmp.refl"

    args = [
        expt_file,
        refl_file,
        "output.experiments=symmetrized.expt",
        "output.reflections=symmetrized.refl",
        "output.html=cosym.html",
        "output.json=cosym.json",
    ]

    if use_known_space_group:
        args.append(f"space_group={space_group.info()}")
    if use_known_lattice_group:
        args.append(f"lattice_group={space_group.info()}")
    if dimensions is not None:
        args.append(f"dimensions={dimensions}")

    dials_cosym.run(args=args)
    assert pathlib.Path("symmetrized.refl").is_file()
    assert pathlib.Path("symmetrized.expt").is_file()
    assert pathlib.Path("cosym.html").is_file()
    assert pathlib.Path("cosym.json").is_file()
    cosym_expts = load.experiment_list("symmetrized.expt", check_format=False)
    assert len(cosym_expts) == len(experiments)
    for expt in cosym_expts:
        if unit_cell is not None:
            assert expt.crystal.get_unit_cell().parameters() == pytest.approx(
                unit_cell.parameters()
            )
        if (
            str(expt.crystal.get_space_group().info()) == "P 6 2 2"
            and str(space_group.info()) == "P 3 2 1"
        ):
            # This is fine
            continue
        assert str(expt.crystal.get_space_group().info()) == str(space_group.info())
        assert expt.crystal.get_space_group() == space_group
