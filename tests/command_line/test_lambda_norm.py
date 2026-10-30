from __future__ import annotations

import json
import shutil
import subprocess

import pytest

from dials.array_family import flex
from dials.command_line import lambda_norm
from dials.util import Sorry


@pytest.fixture
def integrated_nacl(dials_data, run_in_tmp_path):
    processed = dials_data("isis_sxd_nacl_processed")
    image_file = dials_data("isis_sxd_example_data") / "sxd_nacl_run.nxs"
    expt = json.loads((processed / "refined.expt").read_text())
    expt["imageset"][0]["template"] = str(image_file)
    (run_in_tmp_path / "refined.expt").write_text(json.dumps(expt))
    result = subprocess.run(
        [
            shutil.which("dials.tof_integrate"),
            "refined.expt",
            processed / "refined.refl",
        ],
        capture_output=True,
    )
    assert not result.returncode and not result.stderr
    return run_in_tmp_path


def test_lambda_norm(integrated_nacl):
    lambda_norm.run(["integrated.expt", "integrated.refl"])
    for name in ("expt", "refl"):
        assert (integrated_nacl / f"lambda_norm.{name}").is_file()
    for name in ("json", "html"):
        assert not (integrated_nacl / f"lambda_norm.{name}").exists()
    assert (
        "marks the reference wavelength"
        in (integrated_nacl / "dials.lambda_norm.log").read_text()
    )

    before = flex.reflection_table.from_file("integrated.refl")
    after = flex.reflection_table.from_file("lambda_norm.refl")
    w = after["lambda_norm.scale"]
    assert flex.min(w) > 0
    assert (before["intensity.sum.value"] / w).all_approx_equal(
        after["intensity.sum.value"]
    )
    assert (before["intensity.sum.variance"] / (w * w)).all_approx_equal(
        after["intensity.sum.variance"]
    )

    lambda_norm.run(
        [
            "integrated.expt",
            "integrated.refl",
            "output.json=lambda_norm.json",
            "output.html=lambda_norm.html",
        ]
    )
    assert (integrated_nacl / "lambda_norm.html").is_file()
    data = json.loads((integrated_nacl / "lambda_norm.json").read_text())
    assert len(data["training"]["lambda"]) >= 4

    with pytest.raises(Sorry, match="already wavelength-normalised"):
        lambda_norm.run(["lambda_norm.expt", "lambda_norm.refl"])
