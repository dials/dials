from __future__ import annotations

import json
import os
from pathlib import Path

import procrunner

from dials.command_line import stereographic_projection


def test_stereographic_projection(dials_data, tmp_path):
    result = procrunner.run(
        (
            "dials.stereographic_projection",
            dials_data("centroid_test_data", pathlib=True) / "experiments.json",
            "hkl_limit=4",
            "plot.filename=proj.png",
            "json.filename=proj.json",
        ),
        working_directory=tmp_path,
    )
    assert not result.returncode and not result.stderr
    assert tmp_path.joinpath("projections.txt").is_file()
    assert tmp_path.joinpath("proj.png").is_file()
    assert tmp_path.joinpath("proj.json").is_file()

    d = json.loads(tmp_path.joinpath("proj.json").read_text())
    assert set(d) == {"data", "layout"}
    assert d["data"][0]["name"] == "stereographic_projections"
    assert len(d["data"][0]["x"]) == 289


def test_labels(dials_data, tmp_path):
    output_file = tmp_path / "proj.json"
    experiments = sorted(
        dials_data("multi_crystal_proteinase_k", pathlib=True).glob("experiments*.json")
    )
    args = [str(e) for e in experiments] + [
        f"plot.labels={' '.join(str(i) for i in range(len(experiments)))}",
        f"json.filename={output_file}",
        "hkl=1,0,0",
    ]

    cwd = Path.cwd()
    try:
        os.chdir(tmp_path)
        stereographic_projection.run(args)
    finally:
        os.chdir(cwd)
    d = json.loads(output_file.read_bytes())
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
