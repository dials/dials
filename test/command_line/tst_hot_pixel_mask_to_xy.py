from __future__ import division
import os
import libtbx.load_env
from libtbx import easy_run

have_dials_regression = libtbx.env.has_module("dials_regression")
if have_dials_regression:
  dials_regression = libtbx.env.find_in_repositories(
    relative_path="dials_regression",
    test=os.path.isdir)


def exercise_hot_pixel_mask_to_xy():
  if not have_dials_regression:
    print "Skipping exercise_hot_pixel_mask_to_xy(): dials_regression not available."
    return

  data_file = os.path.join(dials_regression, "i23", "hot_pixel_mask",
                           "hot_mask_0.pickle")

  commands = ["dev.dials.hot_pixel_mask_to_xy", data_file]
  command = " ".join(commands)
  result = easy_run.fully_buffered(command=command).raise_if_errors()

  return

def run():
  exercise_hot_pixel_mask_to_xy()
  return

if __name__ == '__main__':
  from dials.test import cd_auto
  with cd_auto(__file__):
    run()
    print "OK"
