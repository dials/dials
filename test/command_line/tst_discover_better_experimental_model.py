from __future__ import division
import os
import libtbx.load_env
from libtbx import easy_run
from libtbx.test_utils import approx_equal
from libtbx.test_utils import open_tmp_directory
from scitbx import matrix

# this import required early to avoid seg fault on some systems
try:
  import scipy.linalg # import dependency
except ImportError, e:
  pass

# apply a random seed to avoid this randomly crashing... I hope
import random
random.seed(12345)

have_dials_regression = libtbx.env.has_module("dials_regression")
if have_dials_regression:
  dials_regression = libtbx.env.find_in_repositories(
    relative_path="dials_regression",
    test=os.path.isdir)

have_xia2_regression = libtbx.env.has_module("xia2_regression")
if have_xia2_regression:
  xia2_regression = libtbx.env.under_build("xia2_regression")


def exercise_1():
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
  assert approx_equal(shift.elems, (-0.242, -0.125, 0.0), eps=1e-2)
  # This produces these two different solutions.
  # The two corresponding origin vectors are:
  # "origin": [ -208.507324496093,   209.20518907699287, -266.11 ]
  # "origin": [ -208.50831812992388, 209.20211805759828, -266.11 ]
  # The remainder of the optimized_datablock.json is identical.
  #
  # TODO: I don't know if both of these are legitimate, or if
  # this is a bug in discover_better_experimental_model.
  os.chdir(cwd)

def exercise_2():

  from dials.test.algorithms.indexing.tst_index import run_one_indexing
  if not have_xia2_regression:
    print "Skipping exercise_2(): xia2_regression not available."
    return

  curdir = os.path.abspath(os.curdir)
  print curdir

  data_dir = os.path.join(xia2_regression, "test_data", "i04_bag_training")

  import glob
  g = glob.glob(os.path.join(data_dir, "*.cbf*"))
  if len(g) == 0:
    print "Skipping exercise_2(): xia2_regression files not downloaded."
    print "Run xia2_regression.fetch_test_data first."
    return

  # beam centre from image headers: 205.28,210.76 mm
  args = ["dials.import", "mosflm_beam_centre=207,212"] + g
  command = " ".join(args)
  #print command
  result = easy_run.fully_buffered(command=command).raise_if_errors()
  assert os.path.exists('datablock.json')

  # spot-finding, just need a subset of the data
  args = ["dials.find_spots", "datablock.json",
          "scan_range=1,10", "scan_range=531,540"]
  command = " ".join(args)
  print command
  result = easy_run.fully_buffered(command=command).raise_if_errors()
  assert os.path.exists('strong.pickle')

  # actually run the beam centre search
  args = ["dials.discover_better_experimental_model", "datablock.json",
          "strong.pickle"]
  command = " ".join(args)
  print command
  result = easy_run.fully_buffered(command=command).raise_if_errors()
  result.show_stdout()
  assert os.path.exists('optimized_datablock.json')

  # look at the results
  from dxtbx.serialize import load
  datablocks = load.datablock("datablock.json", check_format=False)
  original_imageset = datablocks[0].extract_imagesets()[0]
  optimized_datablock = load.datablock('optimized_datablock.json',
                                       check_format=False)
  detector_1 = original_imageset.get_detector()
  detector_2 = optimized_datablock[0].unique_detectors()[0]
  shift = (matrix.col(detector_1[0].get_origin()) -
           matrix.col(detector_2[0].get_origin()))
  print shift

  # check we can actually index the resulting optimized datablock
  from cctbx import uctbx
  expected_unit_cell = uctbx.unit_cell(
    (57.780, 57.800, 150.017, 89.991, 89.990, 90.007))
  expected_rmsds = (0.06, 0.05, 0.001)
  expected_hall_symbol = ' P 1'
  result = run_one_indexing(
    os.path.join(curdir, 'strong.pickle'),
    os.path.join(curdir, 'optimized_datablock.json'), [],
    expected_unit_cell, expected_rmsds, expected_hall_symbol)

def exercise_3():
  if not have_dials_regression:
    print "Skipping exercise(): dials_regression not available."
    return

  data_dir = os.path.join(dials_regression, "indexing_test_data", "phi_scan")
  pickle_path = os.path.join(data_dir, "strong.pickle")
  datablock_path = os.path.join(data_dir, "datablock.json")

  args = ["dials.discover_better_experimental_model",
          datablock_path, pickle_path]

  command = " ".join(args)
  print command
  cwd = os.path.abspath(os.curdir)
  tmp_dir = open_tmp_directory()
  os.chdir(tmp_dir)
  result = easy_run.fully_buffered(command=command).raise_if_errors()
  assert os.path.exists('optimized_datablock.json')
  from dxtbx.serialize import load
  datablocks = load.datablock(datablock_path, check_format=False)
  original_imageset = datablocks[0].extract_imagesets()[0]
  optimized_datablock = load.datablock('optimized_datablock.json',
                                       check_format=False)
  detector_1 = original_imageset.get_detector()
  detector_2 = optimized_datablock[0].unique_detectors()[0]
  shift = (matrix.col(detector_1[0].get_origin()) -
           matrix.col(detector_2[0].get_origin()))
  assert approx_equal(shift.elems, (-1.04661289838, 2.36670867594, 0.0), eps=1e-1)
  os.chdir(cwd)


def run(args):
  exercises = (exercise_1, exercise_2, exercise_3)
  if len(args):
    args = [int(arg) for arg in args]
    for arg in args: assert arg > 0
    exercises = [exercises[arg-1] for arg in args]

  for exercise in exercises:
    exercise()

  print "OK"

if __name__ == '__main__':
  import sys
  from dials.test import cd_auto
  with cd_auto(__file__):
    run(sys.argv[1:])
