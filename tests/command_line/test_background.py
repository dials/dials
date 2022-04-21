from __future__ import annotations

import procrunner


def test(dials_data, tmp_path):
    experiments = dials_data("centroid_test_data", pathlib=True) / "experiments.json"

    result = procrunner.run(
        [
            "dials.background",
            "output.plot=background.png",
            "image=1",
            experiments,
        ],
        working_directory=tmp_path,
    )

    assert not result.returncode and not result.stderr
    assert (tmp_path / "background.png").is_file()

    for line in result.stdout.splitlines():
        if line.startswith(b"Mean background"):
            assert line.endswith(b"0.559")
            break


def test_multiple_imagesets(dials_data, tmp_path):
    filenames = sorted(
        dials_data("thaumatin_grid_scan", pathlib=True).glob("thau_3_2_00*.cbf.bz2")
    )
    filenames.extend(
        sorted(dials_data("centroid_test_data", pathlib=True).glob("centroid_*.cbf"))
    )

    result = procrunner.run(
        [
            "dials.background",
            "output.plot=background.png",
            "image=1,2",
            "size_inches=16,10",
        ]
        + filenames,
        working_directory=tmp_path,
    )

    assert not result.returncode and not result.stderr
    assert (tmp_path / "background.png").is_file()

    lines = result.stdout.splitlines()
    assert b"For imageset 0 image 1:" in lines
    assert b"For imageset 0 image 2:" in lines
    assert b"For imageset 1 image 1:" in lines
    assert b"For imageset 1 image 2:" in lines
