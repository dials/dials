"""Tests for dials.report"""
from __future__ import annotations

import json

import procrunner


def test_report_integrated_data(dials_data, tmpdir):
    """Simple test to check that dials.report completes when given integrated data."""
    data_dir = dials_data("l_cysteine_dials_output", pathlib=True)
    result = procrunner.run(
        [
            "dials.report",
            data_dir / "20_integrated_experiments.json",
            data_dir / "20_integrated.pickle",
        ],
        working_directory=tmpdir,
    )
    assert not result.returncode and not result.stderr
    assert tmpdir.join("dials.report.html").check()


def test_report_scaled_data(dials_data, tmpdir):
    """Test that dials.report works on scaled data."""
    data_dir = dials_data("l_cysteine_4_sweeps_scaled", pathlib=True)
    result = procrunner.run(
        [
            "dials.report",
            data_dir / "scaled_30.refl",
            data_dir / "scaled_30.expt",
            f"json={tmpdir}/report.json",
        ],
        working_directory=tmpdir,
    )
    assert not result.returncode and not result.stderr
    assert tmpdir.join("dials.report.html").check()
    report_json = tmpdir.join("report.json")
    assert report_json.check()
    expected_keys = {
        "strong",
        "centroid",
        "intensity",
        "reference",
        "scan_varying",
        "scaling_model",
        "resolution_graphs",
        "batch_graphs",
        "misc_graphs",
        "scaled_intensity_graphs",
    }
    with report_json.open() as fh:
        d = json.load(fh)
        assert not expected_keys - set(d.keys())
