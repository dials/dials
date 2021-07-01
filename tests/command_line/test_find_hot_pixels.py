import procrunner


def test(dials_data, tmpdir):
    images = dials_data("centroid_test_data").listdir("centroid*.cbf")

    result = procrunner.run(
        [
            "dials.find_spots",
            "nproc=1",
            "output.experiments=spotfinder.expt",
            "output.reflections=spotfinder.refl",
            "output.shoeboxes=True",
        ]
        + images,
        working_directory=tmpdir,
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
        working_directory=tmpdir,
    )
    assert not result.returncode and not result.stderr
    assert tmpdir.join("hot_pixels.mask").check()
    assert (
        b"Found 8 hot pixels" in result.stdout or b"Found 9 hot pixels" in result.stdout
    )
