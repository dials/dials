from __future__ import annotations

import procrunner
import pytest


@pytest.mark.parametrize("dataset", ["centroid_test_data", "thaumatin_grid_scan"])
def test_spot_counts_per_image(dataset, dials_data, tmp_path):
    path = dials_data(dataset, pathlib=True)

    # import the data
    result = procrunner.run(
        ["dials.import", "output.experiments=imported.expt"]
        + list(path.glob("*.cbf*")),
        working_directory=tmp_path,
    )
    assert not result.returncode and not result.stderr
    assert tmp_path.joinpath("imported.expt").is_file()

    # find the spots
    result = procrunner.run(
        ["dials.find_spots", "imported.expt", "nproc=1", "min_spot_size=3"],
        working_directory=tmp_path,
    )
    assert not result.returncode and not result.stderr
    assert tmp_path.joinpath("strong.refl").is_file()

    result = procrunner.run(
        [
            "dials.spot_counts_per_image",
            "imported.expt",
            "strong.refl",
            "plot=spot_counts.png",
        ],
        working_directory=tmp_path,
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
    result = procrunner.run(
        [
            "dials.spot_counts_per_image",
            path / "indexed.expt",
            path / "indexed.refl",
        ],
        working_directory=tmp_path,
    )
    assert result.returncode
    assert b"Only unindexed reflections are currently supported" in result.stderr
