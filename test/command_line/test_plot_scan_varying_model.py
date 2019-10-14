from __future__ import absolute_import, division, print_function

import os
import procrunner


def test(dials_regression, run_in_tmpdir):
    result = procrunner.run(
        [
            "dials.plot_scan_varying_model",
            os.path.join(
                dials_regression,
                "refinement_test_data",
                "multi_sweep_one_sample",
                "glucose_isomerase",
                "SWEEP1",
                "index",
                "sv_refined_experiments.json",
            ),
        ]
    )
    assert not result.returncode and not result.stderr
