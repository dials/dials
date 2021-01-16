from __future__ import absolute_import, division, print_function

import procrunner


def test(dials_data, tmpdir):
    experiments = dials_data("centroid_test_data") / "experiments.json"

    result = procrunner.run(
        [
            "dials.background",
            "output.plot=background.png",
            "image=1",
            str(experiments),
        ],
        working_directory=tmpdir,
    )

    assert not result.returncode and not result.stderr
    assert tmpdir.join("background.png").check()

    for line in result.stdout.splitlines():
        if line.startswith(b"Mean background"):
            assert line.endswith(b"0.426")
            break


def test_file_output_bad_display(dials_data, tmpdir):
    experiments = dials_data("centroid_test_data") / "experiments.json"

    result = procrunner.run(
        [
            "dials.background",
            "output.plot=background.png",
            "image=1",
            str(experiments),
        ],
        working_directory=tmpdir,
        environment_override={"DISPLAY": "invalid_host:0"},
    )

    assert not result.returncode and not result.stderr
    assert tmpdir.join("background.png").check()
