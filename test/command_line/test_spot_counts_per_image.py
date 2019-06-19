from __future__ import absolute_import, division, print_function

import os
from libtbx import easy_run
from glob import glob


def test_spot_counts_per_image(dials_data, run_in_tmpdir):
    path = dials_data("centroid_test_data").strpath

    # import the data
    cmd = "dials.import %s output.experiments=imported.expt" % " ".join(
        glob(os.path.join(path, "*.cbf"))
    )
    easy_run.fully_buffered(cmd).raise_if_errors()
    assert os.path.exists("imported.expt")

    # find the spots
    cmd = "dials.find_spots imported.expt min_spot_size=3"
    easy_run.fully_buffered(cmd).raise_if_errors()
    assert os.path.exists("strong.refl")

    cmd = "dials.spot_counts_per_image imported.expt strong.refl plot=spot_counts.png"
    result = easy_run.fully_buffered(cmd).raise_if_errors()
    assert os.path.exists("spot_counts.png"), result.show_stdout()

    assert (
        "| image | #spots | #spots_no_ice | total_intensity |"
        + " d_min | d_min (distl method 1) | d_min (distl method 2) |"
        in result.stdout_lines
    ), result.stdout_lines
