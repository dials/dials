from __future__ import absolute_import, division, print_function

import procrunner


def test_merge_cbf(dials_data, tmpdir):
    data_dir = dials_data("centroid_test_data")

    g = data_dir.listdir("*.cbf")
    assert len(g) == 9

    cmd = ["dials.merge_cbf", "merge_n_images=3"] + [f.strpath for f in g]
    result = procrunner.run(cmd, working_directory=tmpdir.strpath)
    assert not result.returncode and not result.stderr
    g = tmpdir.listdir("sum_*.cbf", sort=True)
    assert len(g) == 3

    # test alternate mode of accessing image data
    cmd += ["image_prefix=sum2_", "get_raw_data_from_imageset=false"]
    result = procrunner.run(cmd, working_directory=tmpdir.strpath)
    assert not result.returncode and not result.stderr

    g2 = tmpdir.listdir("sum2_*.cbf", sort=True)
    assert len(g2) == 3

    # check summed images are the same in either case
    for f1, f2 in zip(g, g2):
        print("Testing", f1, f2)
        assert f1.read_binary() == f2.read_binary()
