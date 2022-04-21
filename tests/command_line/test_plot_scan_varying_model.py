from __future__ import annotations

import os

import procrunner


def test(dials_regression, run_in_tmp_path):
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
