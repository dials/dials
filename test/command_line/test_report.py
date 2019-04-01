from __future__ import absolute_import, division, print_function

import os
import procrunner


def test_report_integrated_data(dials_regression, run_in_tmpdir):
    """Simple test to check that dials.report completes when given integrated data."""

    result = procrunner.run(
        [
            "dials.report",
            os.path.join(dials_regression, "xia2-28", "20_integrated_experiments.json"),
            os.path.join(dials_regression, "xia2-28", "20_integrated.pickle"),
        ]
    )
    assert result["exitcode"] == 0
    assert result["stderr"] == ""
    assert os.path.exists("dials-report.html")
