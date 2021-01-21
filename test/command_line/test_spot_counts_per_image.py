from __future__ import absolute_import, division, print_function

import procrunner
import pytest


@pytest.mark.parametrize("dataset", ["centroid_test_data", "thaumatin_grid_scan"])
def test_spot_counts_per_image(dataset, dials_data, tmpdir):
    path = dials_data(dataset)

    # import the data
    result = procrunner.run(
        ["dials.import", "output.experiments=imported.expt"] + path.listdir("*.cbf*"),
        working_directory=tmpdir,
    )
    assert not result.returncode and not result.stderr

    assert tmpdir.join("imported.expt").check()

    # find the spots
    result = procrunner.run(
        ["dials.find_spots", "imported.expt", "nproc=1", "min_spot_size=3"],
        working_directory=tmpdir,
    )
    assert not result.returncode and not result.stderr
    assert tmpdir.join("strong.refl").check()

    result = procrunner.run(
        [
            "dials.spot_counts_per_image",
            "imported.expt",
            "strong.refl",
            "plot=spot_counts.png",
        ],
        working_directory=tmpdir,
    )
    assert not result.returncode and not result.stderr
    assert tmpdir.join("spot_counts.png").check(), result.stdout

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
