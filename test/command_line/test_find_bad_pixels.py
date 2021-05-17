import procrunner


def test_find_bad_pixels(dials_data, tmpdir):

    image_files = [f.strpath for f in dials_data("x4wide").listdir("*.cbf", sort=True)]
    result = procrunner.run(
        [
            "dials.find_bad_pixels",
        ]
        + image_files,
        working_directory=tmpdir,
    )
    assert not result.returncode and not result.stderr

    # count bad pixels
    count = 0
    for record in result.stdout.decode().split("\n"):
        if "mask" in record:
            assert record.split()[-1] == "8"
            count += 1

    assert count == 27
