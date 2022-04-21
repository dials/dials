from __future__ import annotations

import procrunner


def test(dials_data, tmp_path):
    input_filename = (
        dials_data("centroid_test_data", pathlib=True) / "imported_experiments.json"
    )

    result = procrunner.run(
        ["dials.estimate_gain", f"input.experiments={input_filename}"],
        working_directory=tmp_path,
    )
    assert not result.returncode and not result.stderr
    assert b"Estimated gain: 1.0" in result.stdout
