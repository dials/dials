from __future__ import absolute_import, division, print_function

import json
import os

import procrunner


def test_stereographic_projection(dials_regression, tmpdir):
    path = os.path.join(dials_regression, "experiment_test_data")
    result = procrunner.run(
        (
            "dials.stereographic_projection",
            "%s/experiment_1.json" % path,
            "hkl_limit=4",
            "plot.show=False",
            "plot.filename=proj.png",
            "json.filename=proj.json",
        ),
        working_directory=tmpdir,
    )
    assert not result.returncode and not result.stderr
    assert tmpdir.join("projections.txt").check()
    assert tmpdir.join("proj.png").check()
    assert tmpdir.join("proj.json").check()

    with tmpdir.join("proj.json").open("rb") as f:
        d = json.load(f)
    assert set(d) == {"data", "layout"}
    assert d["data"][0]["name"] == "stereographic_projections"
    assert len(d["data"][0]["x"]) == 578
