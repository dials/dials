from __future__ import absolute_import, division, print_function

import pytest
import procrunner

from dxtbx.serialize import load


@pytest.mark.parametrize(
    "space_group",
    [
        None,
        pytest.param(
            "P 1", marks=pytest.mark.xfail(reason="cosym bug if setting space_group=P1")
        ),
        "P 4",
    ],
)
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
    experiments = load.experiment_list(
        tmpdir.join("reindexed_experiments.json").strpath, check_format=False
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
