from __future__ import absolute_import, division, print_function

import procrunner


def test(dials_data, tmpdir):
    images = dials_data("centroid_test_data").listdir("centroid*.cbf")

    result = procrunner.run(
        [
            "dials.find_spots",
            "output.experiments=experiments.json",
            "output.reflections=spotfinder.pickle",
            "output.shoeboxes=True",
        ]
        + [f.strpath for f in images],
        working_directory=tmpdir.strpath,
    )
    assert not result["exitcode"] and not result["stderr"]
    assert tmpdir.join("experiments.json").check()
    assert tmpdir.join("spotfinder.pickle").check()

    result = procrunner.run(
        [
            "dials.find_hot_pixels",
            "input.experiments=experiments.json",
            "input.reflections=spotfinder.pickle",
            "output.mask=hot_mask.pickle",
        ],
        working_directory=tmpdir.strpath,
    )
    assert not result["exitcode"] and not result["stderr"]
    assert tmpdir.join("hot_mask.pickle").check()
    assert "Found 8 hot pixels" in result["stdout"]
