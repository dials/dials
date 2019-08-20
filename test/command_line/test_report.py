"""Tests for dials.report"""
from __future__ import absolute_import, division, print_function

import os
import procrunner


def test_report_integrated_data(dials_regression, tmpdir):
    """Simple test to check that dials.report completes when given integrated data."""

    result = procrunner.run(
        [
            "dials.report",
            os.path.join(dials_regression, "xia2-28", "20_integrated_experiments.json"),
            os.path.join(dials_regression, "xia2-28", "20_integrated.pickle"),
        ],
        working_directory=tmpdir,
    )
    assert not result.returncode and not result.stderr
    assert tmpdir.join("dials-report.html").check()


def test_report_scaled_data(dials_data, tmpdir):
    """Test that dials.report works on scaled data."""
    location = dials_data("l_cysteine_4_sweeps_scaled")
    refl = location.join("scaled_30.refl").strpath
    expt = location.join("scaled_30.expt").strpath

    result = procrunner.run(["dials.report", refl, expt], working_directory=tmpdir)
    assert not result.returncode and not result.stderr
    assert tmpdir.join("dials-report.html").check()
