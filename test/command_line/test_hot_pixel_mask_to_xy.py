from __future__ import absolute_import, division, print_function

import os

from libtbx import easy_run

def test_hot_pixel_mask_to_xy(dials_regression, run_in_tmpdir):
  data_file = os.path.join(dials_regression, "i23", "hot_pixel_mask",
                           "hot_mask_0.pickle")
  commands = ["dev.dials.hot_pixel_mask_to_xy", data_file]
  command = " ".join(commands)
  easy_run.fully_buffered(command=command).raise_if_errors()
