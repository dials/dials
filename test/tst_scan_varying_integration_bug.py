from __future__ import division
import glob
import os
from libtbx import easy_run
from libtbx.test_utils import approx_equal, open_tmp_directory
from cctbx import uctbx

import libtbx.load_env
have_xia2_regression = libtbx.env.has_module("xia2_regression")
if have_xia2_regression:
  xia2_regression = libtbx.env.under_build("xia2_regression")
  have_test_data = os.path.exists(os.path.join(xia2_regression, "test_data"))

def exercise_1():
  data_dir = os.path.join(xia2_regression, "test_data", "X4_wide")
  cwd = os.path.abspath(os.curdir)
  tmp_dir = open_tmp_directory()
  os.chdir(tmp_dir)
  print tmp_dir

  g = sorted(glob.glob(os.path.join(data_dir, "X4_wide_M1S4_2_00*.cbf")))
  assert len(g) == 90, len(g)

  commands = [
    "dials.import %s" %" ".join(g),
    "dials.slice_sweep datablock.json scan_range=80,90",
    "dials.find_spots datablock_80_90.json",
    "dials.index datablock_80_90.json strong.pickle space_group=P41212",
    "dials.refine experiments.json indexed.pickle scan_varying=True",
    "dials.integrate refined_experiments.json indexed.pickle",
    "dials.export refined_experiments.json integrated.pickle"
  ]

  for cmd in commands:
    # print cmd
    result = easy_run.fully_buffered(cmd).raise_if_errors()

  integrated_mtz = "hklout.mtz"
  assert os.path.exists(integrated_mtz)
  from iotbx.reflection_file_reader import any_reflection_file
  reader = any_reflection_file(integrated_mtz)
  mtz_object = reader.file_content()
  assert mtz_object.column_labels()[:14] == [
    'H', 'K', 'L', 'M_ISYM', 'BATCH', 'IPR', 'SIGIPR', 'I', 'SIGI',
    'FRACTIONCALC', 'XDET', 'YDET', 'ROT', 'LP']

  assert len(mtz_object.batches()) == 11
  batch = mtz_object.batches()[0]
  expected_unit_cell = uctbx.unit_cell(
    (42.5787, 42.5787, 40.2983, 90, 90, 90))
  assert expected_unit_cell.is_similar_to(uctbx.unit_cell(list(batch.cell()))), (
    expected_unit_cell.parameters(), list(batch.cell()))
  assert mtz_object.space_group().type().lookup_symbol() == "P 41 21 2"
  assert approx_equal(mtz_object.n_reflections(), 7446, eps=2e3)
  os.chdir(cwd)

def run(args):
  if not have_xia2_regression:
    print "Skipping tst_scan_varying_integration_bug.py: xia2_regression not available"
    return

  if not have_test_data:
    print "Skipping tst_scan_varying_integration_bug.py: xia2_regression " + \
          "test data not available. Please run " + \
          "xia2_regression.fetch_test_data first"

    return

  exercises = (exercise_1,)
  if len(args):
    args = [int(arg) for arg in args]
    for arg in args: assert arg > 0
    exercises = [exercises[arg-1] for arg in args]

  for exercise in exercises:
    exercise()

if __name__ == '__main__':
  import sys
  from libtbx.utils import show_times_at_exit
  show_times_at_exit()
  run(sys.argv[1:])
