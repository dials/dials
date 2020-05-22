from __future__ import absolute_import, division, print_function

import json
import os
from dials.command_line import stereographic_projection

import procrunner


def test_stereographic_projection(dials_regression, tmpdir):
    path = os.path.join(dials_regression, "experiment_test_data")
    result = procrunner.run(
        (
            "dials.stereographic_projection",
            "%s/experiment_1.json" % path,
            "hkl_limit=4",
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
    assert len(d["data"][0]["x"]) == 289


def test_labels(dials_data, tmpdir):
    experiments = dials_data("multi_crystal_proteinase_k").listdir(
        fil="experiments*.json", sort=True
    )
    args = [e.strpath for e in experiments] + [
        "plot.labels=%s" % " ".join(str(i) for i in range(len(experiments))),
        "json.filename=proj.json",
        "hkl=1,0,0",
    ]
    with tmpdir.as_cwd():
        stereographic_projection.run(args)
    with tmpdir.join("proj.json").open() as fh:
        d = json.load(fh)
        assert d["data"][0]["hoverinfo"] == "text"
        assert d["data"][0]["hovertext"] == [
            "0",
            "0",
            "1",
            "1",
            "2",
            "2",
            "3",
            "3",
            "4",
            "4",
            "5",
            "5",
            "6",
            "6",
            "7",
            "7",
        ]
