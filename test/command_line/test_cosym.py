from __future__ import absolute_import, division, print_function

import os
import pytest
import procrunner

from cctbx import sgtbx, uctbx
from dxtbx.serialize import load
from dials.array_family import flex
from dials.command_line.cosym import cosym, phil_scope
from dials.algorithms.symmetry.cosym.observers import register_default_cosym_observers
from dials.algorithms.symmetry.cosym._generate_test_data import (
    generate_experiments_reflections,
)
from dials.util.multi_dataset_handling import assign_unique_identifiers


@pytest.mark.parametrize("space_group", [None, "P 1", "P 4"])
def test_cosym(dials_data, tmpdir, space_group):
    mcp = dials_data("multi_crystal_proteinase_k")
    command = ["dials.cosym", "space_group=" + str(space_group)]
    for i in [1, 2, 3, 4, 5, 7, 8, 10]:
        command.append(mcp.join("experiments_%d.json" % i).strpath)
        command.append(mcp.join("reflections_%d.pickle" % i).strpath)

    result = procrunner.run(command, working_directory=tmpdir.strpath)
    assert not result.returncode and not result.stderr
    assert tmpdir.join("symmetrized.refl").check(file=1)
    assert tmpdir.join("symmetrized.expt").check(file=1)
    experiments = load.experiment_list(
        tmpdir.join("symmetrized.expt").strpath, check_format=False
    )
    if space_group is None:
        assert (
            experiments[0].crystal.get_space_group().type().lookup_symbol() == "P 4 2 2"
        )
    else:
        assert (
            experiments[0].crystal.get_space_group().type().lookup_symbol()
            == space_group
        )


def test_cosym_partial_dataset(dials_data, tmpdir):
    """Test how cosym handles partial/bad datasets."""
    mcp = dials_data("multi_crystal_proteinase_k")
    command = ["dials.cosym"]
    for i in [1, 2]:
        command.append(mcp.join("experiments_%d.json" % i).strpath)
        command.append(mcp.join("reflections_%d.pickle" % i).strpath)
    # Make one dataset that will be removed in prefiltering
    r = flex.reflection_table.from_pickle(mcp.join("reflections_8.pickle").strpath)
    r["partiality"] = flex.double(r.size(), 0.1)
    r.as_pickle(tmpdir.join("renamed.refl").strpath)
    command.append(tmpdir.join("renamed.refl").strpath)
    command.append(mcp.join("experiments_8.json").strpath)
    # Add another good dataset at the end of the input list
    command.append(mcp.join("experiments_10.json").strpath)
    command.append(mcp.join("reflections_10.pickle").strpath)

    result = procrunner.run(command, working_directory=tmpdir.strpath)
    assert not result.returncode and not result.stderr
    assert tmpdir.join("symmetrized.refl").check(file=1)
    assert tmpdir.join("symmetrized.expt").check(file=1)
    experiments = load.experiment_list(
        tmpdir.join("symmetrized.expt").strpath, check_format=False
    )
    assert len(experiments) == 3


def test_cosym_partial_dataset_raises_sorry(dials_data, tmpdir):
    """Test how cosym handles partial/bad datasets."""
    mcp = dials_data("multi_crystal_proteinase_k")
    command = ["dials.cosym"]
    command.append(tmpdir.join("renamed.refl").strpath)
    command.append(mcp.join("experiments_8.json").strpath)
    r2 = flex.reflection_table.from_pickle(mcp.join("reflections_10.pickle").strpath)
    r2["partiality"] = flex.double(r2.size(), 0.1)
    r2.as_pickle(tmpdir.join("renamed2.refl").strpath)
    command.append(tmpdir.join("renamed2.refl").strpath)
    command.append(mcp.join("experiments_10.json").strpath)

    result = procrunner.run(command, working_directory=tmpdir.strpath)
    # Sorry exceptions are only raised as text at the system level
    assert result.stderr.startswith(b"Sorry")


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
    tmpdir,
):
    os.chdir(tmpdir.strpath)
    space_group = sgtbx.space_group_info(space_group).group()
    if unit_cell is not None:
        unit_cell = uctbx.unit_cell(unit_cell)
    experiments, reflections, reindexing_ops = generate_experiments_reflections(
        space_group=space_group,
        unit_cell=unit_cell,
        unit_cell_volume=10000,
        sample_size=sample_size,
        map_to_p1=True,
        d_min=1.5,
    )
    params = phil_scope.extract()
    if use_known_space_group:
        params.space_group = space_group.info()
    if use_known_lattice_group:
        params.lattice_group = space_group.info()
    if dimensions is not None:
        params.dimensions = dimensions
    experiments, reflections = assign_unique_identifiers(experiments, reflections)
    cosym_instance = cosym(experiments, reflections, params=params)
    register_default_cosym_observers(cosym_instance)
    cosym_instance.run()
    cosym_instance.export()
    assert os.path.exists(params.output.experiments)
    assert os.path.exists(params.output.reflections)
    assert os.path.exists(params.output.html)
    assert os.path.exists(params.output.json)
    cosym_expts = load.experiment_list(params.output.experiments, check_format=False)
    assert len(cosym_expts) == len(experiments)
    for expt in cosym_expts:
        if unit_cell is not None:
            assert expt.crystal.get_unit_cell().parameters() == pytest.approx(
                unit_cell.parameters()
            )
        assert str(expt.crystal.get_space_group().info()) == str(space_group.info())
        assert expt.crystal.get_space_group() == space_group
