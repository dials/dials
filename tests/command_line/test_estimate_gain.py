import procrunner


def test(dials_data, tmp_path):
    input_filename = dials_data("centroid_test_data").join("datablock.json").strpath

    result = procrunner.run(
        ["dials.estimate_gain", "input.experiments=" + input_filename],
        working_directory=tmp_path,
    )
    assert not result.returncode and not result.stderr
    assert b"Estimated gain: 1.0" in result.stdout
