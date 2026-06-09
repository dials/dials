from __future__ import annotations

import shutil
import subprocess


def test(dials_data, tmp_path):
    input_filename = dials_data("centroid_test_data") / "imported_experiments.json"

    result = subprocess.run(
        [shutil.which("dials.estimate_gain"), f"input.experiments={input_filename}"],
        cwd=tmp_path,
        capture_output=True,
    )
    assert not result.returncode and not result.stderr
    assert b"Estimated gain: 1.0" in result.stdout
