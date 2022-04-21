from __future__ import annotations

import procrunner


def test_merge_cbf(dials_data, tmp_path):
    data_dir = dials_data("centroid_test_data", pathlib=True)

    g = sorted(data_dir.glob("*.cbf"))
    assert len(g) == 9

    cmd = ["dials.merge_cbf", "merge_n_images=3"] + g
    result = procrunner.run(cmd, working_directory=tmp_path)
    assert not result.returncode and not result.stderr
    g = sorted(tmp_path.glob("sum_*.cbf"))
    assert len(g) == 3

    # test alternate mode of accessing image data
    cmd += ["image_prefix=sum2_", "get_raw_data_from_imageset=false"]
    result = procrunner.run(cmd, working_directory=tmp_path)
    assert not result.returncode and not result.stderr

    g2 = sorted(tmp_path.glob("sum2_*.cbf"))
    assert len(g2) == 3

    # check summed images are the same in either case
    for f1, f2 in zip(g, g2):
        print("Testing", f1, f2)
        assert f1.read_bytes() == f2.read_bytes()
