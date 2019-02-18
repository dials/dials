from __future__ import absolute_import, division, print_function

import pytest
import procrunner


@pytest.mark.parametrize("space_group", [None, "P 1"])
def test_cosym(dials_data, tmpdir, space_group):

    mcp = dials_data("multi_crystal_proteinase_k")

    command = ["dials.cosym", "space_group=" + str(space_group)]
    for i in [1, 2, 3, 4, 5, 7, 8, 10]:
        command.append(mcp.join("experiments_%d.json" % i).strpath)
        command.append(mcp.join("reflections_%d.pickle" % i).strpath)

    result = procrunner.run(command, working_directory=tmpdir.strpath)
    assert not result["exitcode"] and not result["stderr"]
    assert tmpdir.join("reindexed_reflections.pickle").check(file=1)
    assert tmpdir.join("reindexed_experiments.json").check(file=1)
