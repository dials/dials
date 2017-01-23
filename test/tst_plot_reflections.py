from __future__ import absolute_import, division
import os
import libtbx.load_env
from libtbx import easy_run

def exercise():
  if not libtbx.env.has_module("dials_regression"):
    print "Skipping test: dials_regression not available"
    return
  dials_regression = libtbx.env.find_in_repositories(
    relative_path="dials_regression",
    test=os.path.isdir)
  data_dir = os.path.join(dials_regression, "indexing_test_data", "i04_weak_data")

  cmd = " ".join(["dials.plot_reflections",
                  os.path.join(data_dir, "datablock_orig.json"),
                  os.path.join(data_dir, "full.pickle"),
                  "scan_range=0,5",
                  ])
  print cmd
  easy_run.fully_buffered(cmd).raise_if_errors()
  assert os.path.exists("centroids.png")

def run():
  exercise()
  print "OK"

if __name__ == '__main__':
  run()
