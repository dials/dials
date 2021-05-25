import procrunner


def test_find_bad_pixels(dials_data, tmp_path):

    image_files = sorted(dials_data("x4wide", pathlib=True).glob("*.cbf"))
    image_files = image_files[:10] + image_files[-10:]
    result = procrunner.run(
        [
            "dials.find_bad_pixels",
        ]
        + image_files,
        working_directory=tmp_path,
    )
    assert not result.returncode and not result.stderr

    # count bad pixels
    count = 0
    for record in result.stdout.decode().split("\n"):
        if "mask" in record:
            assert record.split()[-1] == "8"
            count += 1

    assert count == 23
