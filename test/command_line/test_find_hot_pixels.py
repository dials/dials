from __future__ import absolute_import, division, print_function

import procrunner


def test(dials_data, tmpdir):
    images = dials_data("centroid_test_data").listdir("centroid*.cbf")

    result = procrunner.run(
        [
            "dials.find_spots",
            "output.experiments=spotfinder.expt",
            "output.reflections=spotfinder.refl",
            "output.shoeboxes=True",
        ]
        + [f.strpath for f in images],
        working_directory=tmpdir.strpath,
    )
    assert not result.returncode and not result.stderr
    assert tmpdir.join("spotfinder.expt").check()
    assert tmpdir.join("spotfinder.refl").check()

    result = procrunner.run(
        [
            "dials.find_hot_pixels",
            "input.experiments=spotfinder.expt",
            "input.reflections=spotfinder.refl",
            "output.mask=hot_pixels.mask",
        ],
        working_directory=tmpdir.strpath,
    )
    assert not result.returncode and not result.stderr
    assert tmpdir.join("hot_pixels.mask").check()
    assert (
        "Found 8 hot pixels" in result["stdout"]
        or "Found 9 hot pixels" in result["stdout"]
    )
