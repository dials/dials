"""Tests for dials.report"""
from __future__ import absolute_import, division, print_function

import procrunner


def test_report_integrated_data(dials_data, tmpdir):
    """Simple test to check that dials.report completes when given integrated data."""

    result = procrunner.run(
        [
            "dials.report",
            dials_data("l_cysteine_dials_output") / "20_integrated_experiments.json",
            dials_data("l_cysteine_dials_output") / "20_integrated.pickle",
        ],
        working_directory=tmpdir,
    )
    assert not result.returncode and not result.stderr
    assert tmpdir.join("dials.report.html").check()


def test_report_scaled_data(dials_data, tmpdir):
    """Test that dials.report works on scaled data."""
    result = procrunner.run(
        [
            "dials.report",
            dials_data("l_cysteine_4_sweeps_scaled") / "scaled_30.refl",
            dials_data("l_cysteine_4_sweeps_scaled") / "scaled_30.expt",
        ],
        working_directory=tmpdir,
    )
    assert not result.returncode and not result.stderr
    assert tmpdir.join("dials.report.html").check()
