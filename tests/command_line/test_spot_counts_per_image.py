from __future__ import annotations

import shutil
import subprocess

import pytest


@pytest.mark.parametrize("dataset", ["centroid_test_data", "thaumatin_grid_scan"])
def test_spot_counts_per_image(dataset, dials_data, tmp_path):
    path = dials_data(dataset, pathlib=True)

    # import the data
    result = subprocess.run(
        [shutil.which("dials.import"), "output.experiments=imported.expt"]
        + list(path.glob("*.cbf*")),
        cwd=tmp_path,
        capture_output=True,
    )
    assert not result.returncode and not result.stderr
    assert tmp_path.joinpath("imported.expt").is_file()

    # find the spots
    result = subprocess.run(
        [
            shutil.which("dials.find_spots"),
            "imported.expt",
            "nproc=1",
            "min_spot_size=3",
        ],
        cwd=tmp_path,
        capture_output=True,
    )
    assert not result.returncode and not result.stderr
    assert tmp_path.joinpath("strong.refl").is_file()

    result = subprocess.run(
        [
            shutil.which("dials.spot_counts_per_image"),
            "imported.expt",
            "strong.refl",
            "plot=spot_counts.png",
        ],
        cwd=tmp_path,
        capture_output=True,
    )
    assert not result.returncode and not result.stderr
    assert tmp_path.joinpath("spot_counts.png").is_file()

    assert all(
        s in result.stdout
        for s in (
            b"image",
            b"#spots_no_ice",
            b"total_intensity",
            b"d_min (distl method 1)",
            b"d_min (distl method 2)",
        )
    ), result.stdout


def test_spot_counts_per_image_fails_cleanly_on_indexed(dials_data, tmp_path):
    path = dials_data("insulin_processed", pathlib=True)
    result = subprocess.run(
        [
            shutil.which("dials.spot_counts_per_image"),
            path / "indexed.expt",
            path / "indexed.refl",
        ],
        cwd=tmp_path,
        capture_output=True,
    )
    assert result.returncode
    assert b"Only unindexed reflections are currently supported" in result.stderr
