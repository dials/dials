from __future__ import absolute_import, division, print_function

import glob
import os

import procrunner
import pytest
import scitbx

def test_thing_1(run_in_tmpdir, dials_regression):
  '''Would you like to know more about what this test is supposed to do?
     I would love to. Always remember to use descriptive names.'''

  data_dir = os.path.join(dials_regression, "indexing_test_data", "trypsin")
  pickle_path1 = os.path.join(data_dir, "strong_P1_X6_1_0-1.pickle")
  pickle_path2 = os.path.join(data_dir, "strong_P1_X6_2_0-1.pickle")
  experiments_path1 = os.path.join(data_dir, "datablock_P1_X6_1.json")
  experiments_path2 = os.path.join(data_dir, "datablock_P1_X6_2.json")

  args = ["dials.search_beam_position",
          experiments_path1,
          experiments_path2,
          pickle_path1,
          pickle_path2]

  print(args)
  result = procrunner.run_process(args)
  assert result['stderr'] == '' and result['exitcode'] == 0
  assert os.path.exists('optimized_experiments.json')

  from dxtbx.serialize import load
  experiments = load.experiment_list(experiments_path1, check_format=False)
  original_imageset = experiments.imagesets()[0]
  optimized_experiments = load.experiment_list('optimized_experiments.json',
                                       check_format=False)
  detector_1 = original_imageset.get_detector()
  detector_2 = optimized_experiments.detectors()[0]
  shift = (scitbx.matrix.col(detector_1[0].get_origin()) -
           scitbx.matrix.col(detector_2[0].get_origin()))
  assert shift.elems == pytest.approx((0.037, 0.061, 0.0), abs=1e-1)

def test_thing_2(tmpdir, xia2_regression_build):
  '''Would you like to know more about what this test is supposed to do?
     I would love to. Always remember to use descriptive names.'''

  data_dir = os.path.join(xia2_regression_build, "test_data", "i04_bag_training")

  g = glob.glob(os.path.join(data_dir, "*.cbf*"))
  if not g:
    pytest.skip('xia2_regression files not downloaded. Run xia2_regression.fetch_test_data first')

  tmpdir.chdir()

  # beam centre from image headers: 205.28,210.76 mm
  args = ["dials.import", "mosflm_beam_centre=207,212"] + g
  print(args)
  if os.name != 'nt':
    result = procrunner.run_process(args)
    assert result['stderr'] == '' and result['exitcode'] == 0
  else:
    # Can't run this command on Windows, as it will exceed the maximum Windows command length limits.
    # So, instead:
    import mock
    import sys
    with mock.patch.object(sys, 'argv', args):
      import dials.command_line.dials_import
      dials.command_line.dials_import.Script().run()

  assert os.path.exists('experiments.json')

  # spot-finding, just need a subset of the data
  args = ["dials.find_spots", "experiments.json",
          "scan_range=1,10", "scan_range=531,540"]
  print(args)
  result = procrunner.run_process(args)
  assert result['stderr'] == '' and result['exitcode'] == 0
  assert os.path.exists('strong.pickle')

  # actually run the beam centre search
  args = ["dials.search_beam_position", "experiments.json",
          "strong.pickle"]
  print(args)
  result = procrunner.run_process(args)
  assert result['stderr'] == '' and result['exitcode'] == 0
  assert os.path.exists('optimized_experiments.json')

  # look at the results
  from dxtbx.serialize import load
  experiments = load.experiments("experiments.json", check_format=False)
  original_imageset = experiments.imagesets()[0]
  optimized_experiments = load.experiments('optimized_experiments.json',
                                       check_format=False)
  detector_1 = original_imageset.get_detector()
  detector_2 = optimized_experiments.detectors()[0]
  shift = (scitbx.matrix.col(detector_1[0].get_origin()) -
           scitbx.matrix.col(detector_2[0].get_origin()))
  print(shift)

  # check we can actually index the resulting optimized experiments
  from cctbx import uctbx
  from dials.test.algorithms.indexing.test_index import run_one_indexing
  expected_unit_cell = uctbx.unit_cell(
    (57.780, 57.800, 150.017, 89.991, 89.990, 90.007))
  expected_rmsds = (0.06, 0.05, 0.001)
  expected_hall_symbol = ' P 1'
  run_one_indexing(
    (tmpdir / 'strong.pickle').strpath,
    (tmpdir / 'optimized_experiments.json').strpath, [],
    expected_unit_cell, expected_rmsds, expected_hall_symbol)

def test_thing_3(run_in_tmpdir, dials_regression):
  '''Would you like to know more about what this test is supposed to do?
     I would love to. Always remember to use descriptive names.'''

  data_dir = os.path.join(dials_regression, "indexing_test_data", "phi_scan")
  pickle_path = os.path.join(data_dir, "strong.pickle")
  experiments_path = os.path.join(data_dir, "datablock.json")

  args = ["dials.search_beam_position",
          experiments_path, pickle_path]
  print(args)
  result = procrunner.run_process(args)
  assert result['stderr'] == '' and result['exitcode'] == 0
  assert os.path.exists('optimized_experiments.json')

  from dxtbx.serialize import load
  experiments = load.experiment_list(experiments_path, check_format=False)
  original_imageset = experiments.imagesets()[0]
  optimized_experiments = load.experiment_list('optimized_experiments.json',
                                       check_format=False)
  detector_1 = original_imageset.get_detector()
  detector_2 = optimized_experiments.detectors()[0]
  shift = (scitbx.matrix.col(detector_1[0].get_origin()) -
           scitbx.matrix.col(detector_2[0].get_origin()))
  assert shift.elems == pytest.approx((-0.976, 2.497, 0.0), abs=1e-1)
