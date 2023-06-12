from __future__ import annotations

import shutil
import subprocess


def test(dials_data, tmp_path):
    images = dials_data("centroid_test_data", pathlib=True).glob("centroid*.cbf")

    result = subprocess.run(
        [
            shutil.which("dials.find_spots"),
            "nproc=1",
            "output.experiments=spotfinder.expt",
            "output.reflections=spotfinder.refl",
            "output.shoeboxes=True",
        ]
        + sorted(images),
        cwd=tmp_path,
        capture_output=True,
    )
    assert not result.returncode and not result.stderr
    assert tmp_path.joinpath("spotfinder.expt").is_file()
    assert tmp_path.joinpath("spotfinder.refl").is_file()

    result = subprocess.run(
        [
            "dials.find_hot_pixels",
            "input.experiments=spotfinder.expt",
            "input.reflections=spotfinder.refl",
            "output.mask=hot_pixels.mask",
        ],
        cwd=tmp_path,
        capture_output=True,
    )
    assert not result.returncode and not result.stderr
    assert tmp_path.joinpath("hot_pixels.mask").is_file()
    assert (
        b"Found 8 hot pixels" in result.stdout or b"Found 9 hot pixels" in result.stdout
    )
