from __future__ import division
import os
import libtbx.load_env
from libtbx import easy_run
from libtbx.test_utils import approx_equal
from libtbx.test_utils import open_tmp_directory
from scitbx import matrix

# apply a random seed to avoid this randomly crashing... I hope
import random
random.seed(12345)

have_dials_regression = libtbx.env.has_module("dials_regression")
if have_dials_regression:
  dials_regression = libtbx.env.find_in_repositories(
    relative_path="dials_regression",
    test=os.path.isdir)

def exercise():
  if not have_dials_regression:
    print "Skipping exercise(): dials_regression not available."
    return

  data_dir = os.path.join(dials_regression, "indexing_test_data", "trypsin")
  pickle_path1 = os.path.join(data_dir, "strong_P1_X6_1_0-1.pickle")
  pickle_path2 = os.path.join(data_dir, "strong_P1_X6_2_0-1.pickle")
  datablock_path1 = os.path.join(data_dir, "datablock_P1_X6_1.json")
  datablock_path2 = os.path.join(data_dir, "datablock_P1_X6_2.json")

  args = ["dials.discover_better_experimental_model",
          datablock_path1,
          datablock_path2,
          pickle_path1,
          pickle_path2]

  command = " ".join(args)
  print command
  cwd = os.path.abspath(os.curdir)
  tmp_dir = open_tmp_directory()
  os.chdir(tmp_dir)
  result = easy_run.fully_buffered(command=command).raise_if_errors()
  assert os.path.exists('optimized_datablock.json')
  from dxtbx.serialize import load
  datablocks = load.datablock(datablock_path1, check_format=False)
  original_imageset = datablocks[0].extract_imagesets()[0]
  optimized_datablock = load.datablock('optimized_datablock.json',
                                       check_format=False)
  detector_1 = original_imageset.get_detector()
  detector_2 = optimized_datablock[0].unique_detectors()[0]
  shift = (matrix.col(detector_1[0].get_origin()) -
           matrix.col(detector_2[0].get_origin()))
  print shift.elems
  print (-0.178, -0.041, 0.0)
  assert approx_equal(shift.elems, (-0.178, -0.041, 0.0), eps=1e-2)
  os.chdir(cwd)

def run():
  exercise()
  print "OK"

if __name__ == '__main__':
  from dials.test import cd_auto
  with cd_auto(__file__):
    run()
