from __future__ import annotations

import procrunner


def test(dials_data, run_in_tmp_path):
    experiments = (
        dials_data("refinement_test_data", pathlib=True)
        / "glucose_isomerase_sv_refined.json"
    )
    result = procrunner.run(
        [
            "dials.plot_scan_varying_model",
            experiments,
        ]
    )
    assert not result.returncode and not result.stderr
