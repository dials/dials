from __future__ import annotations

import json
import shutil
import subprocess

import pytest

from dials.array_family import flex


@pytest.fixture(scope="module")
def tof_integration_input(dials_data, tmp_path_factory):
    tmp_path = tmp_path_factory.mktemp("tof_integration_input")
    image_file = dials_data("isis_sxd_example_data") / "sxd_nacl_run.nxs"
    processed = dials_data("isis_sxd_nacl_processed")

    # Update .expt file with correct image path
    with open(processed / "refined.expt", "rb") as g:
        expt_json = json.load(g)
    expt_json["imageset"][0]["template"] = str(image_file)
    expt_path = tmp_path / "refined.expt"
    with open(expt_path, "w") as g:
        json.dump(expt_json, g)

    # Reduce number of reflections
    reflections = flex.reflection_table.from_msgpack_file(processed / "refined.refl")
    reflections = reflections.select(reflections["intensity.sum.value"] > 30000)
    refl_path = tmp_path / "refined.refl"
    reflections.as_msgpack_file(refl_path)

    return expt_path, refl_path


# Each case is a separate test so that they can run in parallel under xdist
@pytest.mark.parametrize(
    "option",
    [
        pytest.param(None, id="summation"),
        pytest.param("mask=seed_skewness", id="summation_seed_skewness"),
        "method=profile_1d_ibix",
        "method=profile_1d_ic",  # Ikeda-Carpenter
        "method=profile_3d_gutmann",
        "method=profile_3d_ic",  # Ikeda-Carpenter x bivariate Gaussian
        "method=profile_3d_ibix",  # iBIX x bivariate Gaussian
    ],
)
def test_tof_integration(tof_integration_input, tmp_path, option):
    expt_path, refl_path = tof_integration_input
    command = [shutil.which("dials.tof_integrate"), expt_path, refl_path]
    if option:
        command.append(option)
    result = subprocess.run(command, cwd=tmp_path, capture_output=True)
    assert not result.returncode and not result.stderr
