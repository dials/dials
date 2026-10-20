"""Tests for dials.report"""

from __future__ import annotations

import json
import shutil
import subprocess


def test_report_integrated_data(dials_data, tmp_path):
    """Simple test to check that dials.report completes when given integrated data."""
    data_dir = dials_data("l_cysteine_dials_output")
    result = subprocess.run(
        [
            shutil.which("dials.report"),
            data_dir / "20_integrated_experiments.json",
            data_dir / "20_integrated.pickle",
        ],
        cwd=tmp_path,
        capture_output=True,
    )
    assert not result.returncode and not result.stderr
    assert (tmp_path / "dials.report.html").is_file()


def test_report_scaled_data(dials_data, tmp_path):
    """Test that dials.report works on scaled data."""
    data_dir = dials_data("l_cysteine_4_sweeps_scaled")
    result = subprocess.run(
        [
            shutil.which("dials.report"),
            data_dir / "scaled_30.refl",
            data_dir / "scaled_30.expt",
            f"json={tmp_path}/report.json",
        ],
        cwd=tmp_path,
        capture_output=True,
    )
    assert not result.returncode and not result.stderr
    assert (tmp_path / "dials.report.html").is_file()
    report_json = tmp_path / "report.json"
    assert report_json.is_file()
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
        "scaling_tables",
        "xtriage_output",
        "image_range_tables",
    }
    with report_json.open(encoding="utf-8") as fh:
        d = json.load(fh)
        assert not expected_keys - set(d.keys())
