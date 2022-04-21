from __future__ import annotations

import procrunner


def test(dials_data, tmp_path):
    images = dials_data("centroid_test_data", pathlib=True).glob("centroid*.cbf")

    result = procrunner.run(
        [
            "dials.find_spots",
            "nproc=1",
            "output.experiments=spotfinder.expt",
            "output.reflections=spotfinder.refl",
            "output.shoeboxes=True",
        ]
        + sorted(images),
        working_directory=tmp_path,
    )
    assert not result.returncode and not result.stderr
    assert tmp_path.joinpath("spotfinder.expt").is_file()
    assert tmp_path.joinpath("spotfinder.refl").is_file()

    result = procrunner.run(
        [
            "dials.find_hot_pixels",
            "input.experiments=spotfinder.expt",
            "input.reflections=spotfinder.refl",
            "output.mask=hot_pixels.mask",
        ],
        working_directory=tmp_path,
    )
    assert not result.returncode and not result.stderr
    assert tmp_path.joinpath("hot_pixels.mask").is_file()
    assert (
        b"Found 8 hot pixels" in result.stdout or b"Found 9 hot pixels" in result.stdout
    )
