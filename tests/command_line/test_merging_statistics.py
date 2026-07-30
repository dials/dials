from __future__ import annotations

import pathlib
import shutil
import subprocess


def test_merging_statistics(dials_data, tmp_path):
    """Test the command line script with LCY data"""
    location = dials_data("l_cysteine_4_sweeps_scaled")
    refls = location / "scaled_20_25.refl"
    expts = location / "scaled_20_25.expt"

    command = [
        shutil.which("dials.merging_statistics"),
        refls,
        expts,
        "calculate_sigma_tau_stats=True",
        "additional_stats=True",
    ]
    result = subprocess.run(command, cwd=tmp_path, capture_output=True)
    assert not result.returncode and not result.stderr
    assert (tmp_path / "dials.merging_statistics.html").is_file()
    merge_json = tmp_path / "dials.merging_statistics.json"
    assert merge_json.is_file()


def test_merging_statistics_mtz(dials_data, tmp_path):
    """Test the command line script with LCY data"""
    location = dials_data("l_cysteine_4_sweeps_scaled")
    refls = location / "scaled_20_25.refl"
    expts = location / "scaled_20_25.expt"

    # First export the data
    command = [shutil.which("dials.export"), refls, expts]
    result = subprocess.run(command, cwd=tmp_path, capture_output=True)
    assert not result.returncode and not result.stderr
    assert pathlib.Path(tmp_path / "scaled.mtz").is_file()

    command = [
        shutil.which("dials.merging_statistics"),
        "scaled.mtz",
        "calculate_sigma_tau_stats=True",
        "additional_stats=True",
    ]
    result = subprocess.run(command, cwd=tmp_path, capture_output=True)
    assert not result.returncode and not result.stderr
    assert (tmp_path / "dials.merging_statistics.html").is_file()
    merge_json = tmp_path / "dials.merging_statistics.json"
    assert merge_json.is_file()
