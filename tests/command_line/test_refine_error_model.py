from __future__ import annotations

import json
import shutil
import subprocess

import pytest


@pytest.mark.parametrize("grouping", ["combined", "individual", "grouped"])
def test_standalone_error_model_refinement_on_scaled_data(
    dials_data, tmp_path, grouping
):
    location = dials_data("l_cysteine_4_sweeps_scaled", pathlib=True)
    refls = location / "scaled_20_25.refl"
    expts = location / "scaled_20_25.expt"

    # Use a range of options, some default.
    command = [
        shutil.which("dials.refine_error_model"),
        refls,
        expts,
        "json=error_model.json",
        "html=error_model.html",
        "intensity_choice=combine",
        "combine.Imid=250",
        "basic.minimisation=individual",
        f"grouping={grouping}",
    ]
    if grouping == "grouped":
        command.append("error_model_group='0'")
        command.append("error_model_group='1'")

    result = subprocess.run(command, cwd=tmp_path, capture_output=True)
    assert not result.returncode and not result.stderr
    assert tmp_path.joinpath("error_model.html").is_file()
    assert tmp_path.joinpath("error_model.json").is_file()
    with open(tmp_path.joinpath("error_model.json"), "r") as f:
        d = json.load(f)
    n_models = sum(
        "normal_distribution_plot" in k for k in d["error_model_plots"].keys()
    )
    if grouping == "combined":
        assert n_models == 1
    else:
        assert n_models == 2
