from __future__ import annotations

import shutil
import subprocess


def test_run(dials_data, tmp_path):
    subprocess.run(
        (
            shutil.which("dials.plot_reflections"),
            dials_data("centroid_test_data", pathlib=True) / "experiments.json",
            dials_data("centroid_test_data", pathlib=True) / "integrated.refl",
            "scan_range=0,5",
        ),
        cwd=tmp_path,
        capture_output=True,
    ).check_returncode()
    assert (tmp_path / "centroids.png").is_file()
