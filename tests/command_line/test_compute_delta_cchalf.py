"""Tests for dials.compute_delta_cchalf."""

from __future__ import annotations

import os

import procrunner
import pytest

from dials.command_line.compute_delta_cchalf import CCHalfFromMTZ, phil_scope, run


def check_cchalf_result(fileobj):
    """Inspect the result"""
    lines = fileobj.readlines()
    assert lines[0] == "1 -0.004673\n"
    assert lines[1] == "0 0.001234\n"


def test_suitable_exit_on_bad_input(dials_data, run_in_tmp_path):
    location = dials_data("l_cysteine_4_sweeps_scaled", pathlib=True)
    refl = location / "scaled_35.refl"
    expt = location / "scaled_35.expt"

    args = [str(refl), str(expt)]
    with pytest.raises(SystemExit):
        run(args)

    args = [str(refl), str(expt), "mode=image_group", "group_size=10000"]
    with pytest.raises(SystemExit):
        run(args)


def test_compute_delta_cchalf_scaled_data(dials_data, tmp_path):
    """Test dials.compute_delta_cchalf on scaled data."""
    location = dials_data("l_cysteine_4_sweeps_scaled", pathlib=True)
    refls = location / "scaled_20_25.refl"
    expts = location / "scaled_20_25.expt"

    # set cutoff to 0.0 to force one to be 'rejected'
    command = [
        "dials.compute_delta_cchalf",
        refls,
        expts,
        "stdcutoff=0.0",
        "output.reflections=filtered.refl",
        "output.experiments=filtered.expt",
    ]
    result = procrunner.run(command, working_directory=tmp_path)
    assert not result.returncode and not result.stderr
    assert (tmp_path / "filtered.expt").is_file()
    assert (tmp_path / "filtered.refl").is_file()
    assert (tmp_path / "delta_cchalf.dat").is_file()
    assert (tmp_path / "compute_delta_cchalf.html").is_file()
    with open(tmp_path / "delta_cchalf.dat") as f:
        check_cchalf_result(f)


def test_compute_delta_cchalf_scaled_data_mtz(dials_data, tmp_path):
    """Test dials.compute_delta_cchalf on scaled data."""
    location = dials_data("l_cysteine_4_sweeps_scaled", pathlib=True)
    refls = location / "scaled_20_25.refl"
    expts = location / "scaled_20_25.expt"

    # First export the data
    command = ["dials.export", refls, expts, "partiality_threshold=0.99"]
    result = procrunner.run(command, working_directory=tmp_path)
    assert not result.returncode and not result.stderr
    assert (tmp_path / "scaled.mtz").is_file()

    # set cutoff to 0.0 to force one to be 'rejected'
    command = [
        "dials.compute_delta_cchalf",
        f"mtzfile={tmp_path / 'scaled.mtz'}",
        "stdcutoff=0.0",
    ]
    result = procrunner.run(command, working_directory=tmp_path)
    assert not result.returncode and not result.stderr
    assert (tmp_path / "delta_cchalf.dat").is_file()
    assert (tmp_path / "compute_delta_cchalf.html").is_file()
    with open(tmp_path / "delta_cchalf.dat") as f:
        check_cchalf_result(f)


def test_compute_delta_cchalf(dials_regression):
    """Test compute delta cchalf on an integrated mtz."""

    filename = os.path.join(
        dials_regression, "delta_cchalf_test_data", "test.XDS_ASCII.mtz"
    )
    params = phil_scope.extract()
    params.nbins = 1

    script = CCHalfFromMTZ(params, filename)
    script.run()

    cchalf_i = script.results_summary["per_dataset_delta_cc_half_values"][
        "delta_cc_half_values"
    ]

    assert abs(100 * script.results_summary["mean_cc_half"] - 94.582) < 1e-3
    assert (
        abs(100 * script.results_summary["mean_cc_half"] - 100 * cchalf_i[1] - 79.587)
        < 1e-3
    )
    assert (
        abs(100 * script.results_summary["mean_cc_half"] - 100 * cchalf_i[0] - 94.238)
        < 1e-3
    )
