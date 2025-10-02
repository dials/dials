from __future__ import annotations

import shutil
import subprocess

import pytest

from dials.array_family import flex


def test_against_dials_integrate(dials_data, tmp_path):
    ## ensure insulin folder exists
    dials_data("insulin", pathlib=True)

    # Run as a single job to avoid splitting reflections
    subprocess.run(
        (
            shutil.which("dials.integrate"),
            dials_data("insulin_processed", pathlib=True) / "refined.expt",
            dials_data("insulin_processed", pathlib=True) / "refined.refl",
            "mp.njobs=1",
            "mp.nproc=1",
            "scan_range=1,2",
        ),
        cwd=tmp_path,
        capture_output=True,
    ).check_returncode()

    subprocess.run(
        (
            shutil.which("dev.dials.simple_integrate"),
            dials_data("insulin_processed", pathlib=True) / "refined.expt",
            dials_data("insulin_processed", pathlib=True) / "refined.refl",
            "scan_range=1,2",
        ),
        cwd=tmp_path,
        capture_output=True,
    ).check_returncode()

    simple_refl = flex.reflection_table.from_file(tmp_path / "simple_integrated.refl")
    dials_refl = flex.reflection_table.from_file(tmp_path / "integrated.refl")

    matches = dials_refl.match_with_reference(simple_refl)[0]
    dials_refl = dials_refl.select(matches)

    assert len(simple_refl) == len(dials_refl)
    assert sum(
        simple_refl["intensity.sum.value"] - dials_refl["intensity.sum.value"]
    ) == pytest.approx(0.0, abs=1e-7)
    assert sum(
        simple_refl["intensity.prf.value"] - dials_refl["intensity.prf.value"]
    ) == pytest.approx(0.0, abs=1e-7)
    assert sum(
        simple_refl["background.sum.value"] - dials_refl["background.sum.value"]
    ) == pytest.approx(0.0, abs=1e-7)
    assert sum(
        simple_refl["background.mean"] - dials_refl["background.mean"]
    ) == pytest.approx(0.0, abs=1e-7)
    assert sum(
        simple_refl["background.dispersion"] - dials_refl["background.dispersion"]
    ) == pytest.approx(0.0, abs=1e-7)
    assert sum(
        simple_refl["background.mse"] - dials_refl["background.mse"]
    ) == pytest.approx(0.0, abs=1e-7)
    assert sum(
        simple_refl["background.sum.variance"] - dials_refl["background.sum.variance"]
    ) == pytest.approx(0.0, abs=1e-7)
    assert sum(
        simple_refl["intensity.prf.variance"] - dials_refl["intensity.prf.variance"]
    ) == pytest.approx(0.0, abs=1e-7)
    assert sum(
        simple_refl["intensity.sum.variance"] - dials_refl["intensity.sum.variance"]
    ) == pytest.approx(0.0, abs=1e-7)
    assert sum(simple_refl["lp"] - dials_refl["lp"]) == pytest.approx(0.0, abs=1e-7)
    assert sum(simple_refl["zeta"] - dials_refl["zeta"]) == pytest.approx(0.0, abs=1e-7)
    assert sum(simple_refl["d"] - dials_refl["d"]) == pytest.approx(0.0, abs=1e-7)
    assert sum(simple_refl["partiality"] - dials_refl["partiality"]) == pytest.approx(
        0.0, abs=1e-7
    )
    assert (
        simple_refl["num_pixels.foreground"] == dials_refl["num_pixels.foreground"]
    ).count(True) == len(simple_refl)
    assert (
        simple_refl["num_pixels.background"] == dials_refl["num_pixels.background"]
    ).count(True) == len(simple_refl)
    assert (simple_refl["num_pixels.valid"] == dials_refl["num_pixels.valid"]).count(
        True
    ) == len(simple_refl)
    assert (
        simple_refl["num_pixels.background_used"]
        == dials_refl["num_pixels.background_used"]
    ).count(True) == len(simple_refl)
    assert sum(
        simple_refl["profile.correlation"] - dials_refl["profile.correlation"]
    ) == pytest.approx(0.0, abs=1e-7)
