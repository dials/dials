from __future__ import absolute_import, division, print_function

import pytest
import procrunner

from dxtbx.serialize import load
from dials.array_family import flex


@pytest.mark.parametrize("space_group", [None, "P 1", "P 4"])
def test_cosym(dials_data, tmpdir, space_group):
    mcp = dials_data("multi_crystal_proteinase_k")
    command = ["dials.cosym", "space_group=" + str(space_group)]
    for i in [1, 2, 3, 4, 5, 7, 8, 10]:
        command.append(mcp.join("experiments_%d.json" % i).strpath)
        command.append(mcp.join("reflections_%d.pickle" % i).strpath)

    result = procrunner.run(command, working_directory=tmpdir.strpath)
    assert not result["exitcode"] and not result["stderr"]
    assert tmpdir.join("reindexed_reflections.refl").check(file=1)
    assert tmpdir.join("reindexed_experiments.expt").check(file=1)
    experiments = load.experiment_list(
        tmpdir.join("reindexed_experiments.expt").strpath, check_format=False
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
    assert not result["exitcode"] and not result["stderr"]
    assert tmpdir.join("reindexed_reflections.refl").check(file=1)
    assert tmpdir.join("reindexed_experiments.expt").check(file=1)
    experiments = load.experiment_list(
        tmpdir.join("reindexed_experiments.expt").strpath, check_format=False
    )
    assert len(experiments) == 3

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
    assert result["stderr"][0:5] == "Sorry"
