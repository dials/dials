from __future__ import annotations

import os
import shutil
import subprocess


def test(dials_data, tmp_path):
    experiments = (
        dials_data("refinement_test_data", pathlib=True)
        / "glucose_isomerase_sv_refined.json"
    )
    env = os.environ.copy()
    env[
        "PYTHONDEVMODE"
    ] = ""  # Temporarily disable a developermode warning from pyparsing from mathtext in matplotlib. Try removing after June 2023
    result = subprocess.run(
        [
            shutil.which("dials.plot_scan_varying_model"),
            experiments,
        ],
        cwd=tmp_path,
        env=env,
        capture_output=True,
    )
    assert not result.returncode and not result.stderr
