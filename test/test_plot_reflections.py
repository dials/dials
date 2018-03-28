from __future__ import absolute_import, division, print_function

import os

from libtbx import easy_run

def test_run(dials_regression, tmpdir):
  tmpdir.chdir()
  data_dir = os.path.join(dials_regression, "indexing_test_data", "i04_weak_data")

  cmd = " ".join(["dials.plot_reflections",
                  os.path.join(data_dir, "datablock_orig.json"),
                  os.path.join(data_dir, "full.pickle"),
                  "scan_range=0,5",
                  ])
  print(cmd)
  easy_run.fully_buffered(cmd).raise_if_errors()
  assert os.path.exists("centroids.png")
