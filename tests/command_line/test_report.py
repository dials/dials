"""Tests for dials.report"""

from __future__ import annotations

import json
import re
import shutil
import subprocess

from dxtbx.model.experiment_list import ExperimentListFactory


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


def test_report_tof_refined_data(dials_data, tmp_path):
    """Test that dials.report works on time-of-flight refined data, with the scan
    properties shown in place of oscillation."""
    data_dir = dials_data("isis_sxd_nacl_processed")
    result = subprocess.run(
        [
            shutil.which("dials.report"),
            data_dir / "refined.expt",
            data_dir / "refined.refl",
            f"json={tmp_path}/report.json",
        ],
        cwd=tmp_path,
        capture_output=True,
    )
    assert not result.returncode and not result.stderr
    report_html = tmp_path / "dials.report.html"
    assert report_html.is_file()

    # Scan table shows the image range and the range of each scan property
    scan = ExperimentListFactory.from_json_file(
        data_dir / "refined.expt", check_format=False
    )[0].scan
    tof = scan.get_property("time_of_flight")
    tof_bin_widths = scan.get_property("time_of_flight_bin_widths")
    html = report_html.read_text(encoding="utf-8")
    cells = [c.strip() for c in re.findall(r"<td>(.*?)</td>", html, re.DOTALL)]
    scan_cells = cells[cells.index("<strong>Scan:</strong>") :]
    assert scan_cells[1:5] == [
        "Image range:",
        "%i, %i" % scan.get_image_range(),
        "Time of flight:",
        f"{tof[0]:.3f} - {tof[-1]:.3f}",
    ]
    assert scan_cells[5:8] == [
        "",
        "Time of flight bin widths:",
        f"{tof_bin_widths[0]:.3f} - {tof_bin_widths[-1]:.3f}",
    ]
    assert "Oscillation:" not in cells

    # Centroid differences are shown vs time of flight rather than phi
    with (tmp_path / "report.json").open(encoding="utf-8") as fh:
        centroid = json.load(fh)["centroid"]
    assert {"centroid_mean_differences_vs_tof", "centroid_rmsd_vs_tof"} <= set(centroid)
    assert not {
        "centroid_mean_differences_vs_phi",
        "centroid_rmsd_vs_phi",
    } & set(centroid)
