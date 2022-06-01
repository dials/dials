from __future__ import annotations

import pickle

import procrunner


def return_locations():
    locations = {
        (1042, 980),
        (1042, 988),
        (1042, 1474),
        (1060, 980),
        (1060, 988),
        (1060, 1474),
        (1060, 1482),
        (1253, 980),
        (1254, 980),
        (1254, 988),
        (1254, 1474),
        (1254, 1482),
        (1272, 1474),
        (1272, 1482),
        (1285, 1235),
        (1285, 1237),
        (1465, 980),
        (1466, 980),
        (1466, 988),
        (1466, 1474),
        (1484, 988),
        (1484, 1482),
        (1696, 1482),
    }
    return locations


def test_find_bad_pixels(dials_data, tmp_path):

    image_files = sorted(dials_data("x4wide", pathlib=True).glob("*.cbf"))
    image_files = image_files[:10] + image_files[-10:]
    result = procrunner.run(
        [
            "dials.find_bad_pixels",
            "mask=pixels.mask",
        ]
        + image_files,
        working_directory=tmp_path,
    )
    assert not result.returncode and not result.stderr

    locations = return_locations()
    # verify bad pixels
    for record in result.stdout.decode().split("\n"):
        if "mask" in record and "Wrote" not in record:
            tokens = tuple(map(int, record.split()[1:]))
            assert tokens[-1] == 16
            position = tokens[:2]
            locations.remove(position)

    assert len(locations) == 0

    assert (tmp_path / "pixels.mask").is_file()
    with open(tmp_path / "pixels.mask", "rb") as fh:
        mask = pickle.load(fh, encoding="bytes")

    # X4-wide image size = 2463, 2527
    n_fast = 2463
    locations = return_locations()
    # locations are reported slow fast

    assert mask.count(False) == len(locations)
    for loc in locations:
        idx = (loc[0] * n_fast) + loc[1]
        assert mask[idx] is False

    # check that this mask can be used in import
    image_files = sorted(dials_data("x4wide", pathlib=True).glob("*.cbf"))[:3]
    result = procrunner.run(
        [
            "dials.import",
            "mask=pixels.mask",
        ]
        + image_files,
        working_directory=tmp_path,
    )
    assert not result.returncode and not result.stderr
