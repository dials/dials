from __future__ import absolute_import, division, print_function

import procrunner


def test(dials_data, tmpdir):
    images = dials_data("centroid_test_data").listdir("centroid*.cbf")

    result = procrunner.run(
        [
            "dials.find_spots",
            "output.experiments=experiments.expt",
            "output.reflections=spotfinder.refl",
            "output.shoeboxes=True",
        ]
        + [f.strpath for f in images],
        working_directory=tmpdir.strpath,
    )
    assert not result["exitcode"] and not result["stderr"]
    assert tmpdir.join("experiments.expt").check()
    assert tmpdir.join("spotfinder.refl").check()

    result = procrunner.run(
        [
            "dials.find_hot_pixels",
            "input.experiments=experiments.expt",
            "input.reflections=spotfinder.refl",
            "output.mask=hot_mask.pickle",
        ],
        working_directory=tmpdir.strpath,
    )
    assert not result["exitcode"] and not result["stderr"]
    assert tmpdir.join("hot_mask.pickle").check()
    assert (
        "Found 8 hot pixels" in result["stdout"]
        or "Found 9 hot pixels" in result["stdout"]
    )
