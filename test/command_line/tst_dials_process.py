from __future__ import division
import glob
import os
from libtbx import easy_run
from libtbx.test_utils import approx_equal, open_tmp_directory
from cctbx import uctbx

import libtbx.load_env
have_dials_regression = libtbx.env.has_module("dials_regression")
if have_dials_regression:
  dials_regression = libtbx.env.find_in_repositories(
    relative_path="dials_regression",
    test=os.path.isdir)

def exercise_1():
  data_dir = os.path.join(dials_regression, "mosflm_hg_images")
  cwd = os.path.abspath(os.curdir)
  tmp_dir = open_tmp_directory()
  os.chdir(tmp_dir)
  print tmp_dir
  g = sorted(glob.glob(os.path.join(data_dir, "hg_*.mar1600")))
  assert len(g) == 84
  hall_symbol =  '-R 3 2"'
  cmd = " ".join(["dev.dials.process",
                  "refinement.parameterisation.crystal.scan_varying=True",
                  "indexing.known_symmetry.space_group=\"Hall: %s\"" %hall_symbol.replace('"', '\\"'),
                  "maximum_spot_error=3",
                  "maximum_phi_error=2",
                  "integration.block.threshold=0.99",
                  "template=%s" % os.path.join(data_dir, "hg_###.mar1600"),
                  ]
                 )
  # print cmd
  result = easy_run.fully_buffered(cmd).raise_if_errors()
  for out_file in ['datablock.json', 'refined_experiments.json',
                   'integrated.mtz', 'integrated.pickle', 'strong.pickle']:
    assert os.path.exists(out_file)

  from iotbx.reflection_file_reader import any_reflection_file
  reader = any_reflection_file('integrated.mtz')
  mtz_object = reader.file_content()
  assert mtz_object.column_labels() == [
    'H', 'K', 'L', 'M_ISYM', 'BATCH', 'IPR', 'SIGIPR', 'I', 'SIGI',
    'FRACTIONCALC', 'XDET', 'YDET', 'ROT', 'LP', 'DQE']

  assert len(mtz_object.batches()) == 84
  batch = mtz_object.batches()[0]
  expected_unit_cell = uctbx.unit_cell(
    (58.373, 58.373, 155.939, 90, 90, 120))
  assert expected_unit_cell.is_similar_to(uctbx.unit_cell(list(batch.cell())))
  assert mtz_object.space_group().type().hall_symbol() == hall_symbol
  assert approx_equal(mtz_object.n_reflections(), 24465, eps=2e3)
  os.chdir(cwd)


def exercise_2():
  data_dir = os.path.join(dials_regression, "xia2_demo_data")
  cwd = os.path.abspath(os.curdir)
  tmp_dir = open_tmp_directory()
  os.chdir(tmp_dir)
  g = sorted(glob.glob(os.path.join(data_dir, "insulin*.img")))
  assert len(g) == 45
  hall_symbol =  " I 2 2 3"
  cmd = " ".join(["dev.dials.process",
                  "refinement.parameterisation.crystal.scan_varying=True",
                  "indexing.known_symmetry.space_group=\"Hall: %s\"" %hall_symbol,
                  "maximum_spot_error=3",
                  "maximum_phi_error=2",
                  "integration.block.threshold=0.99",
                  "template=%s" % os.path.join(data_dir, "insulin_1_###.img"),
                  ]
                + ['"%s"' %p for p in g]
                 )
  #print cmd
  result = easy_run.fully_buffered(cmd).raise_if_errors()
  for out_file in ['datablock.json', 'refined_experiments.json',
                   'integrated.mtz', 'integrated.pickle', 'strong.pickle']:
    assert os.path.exists(out_file)

  from iotbx.reflection_file_reader import any_reflection_file
  reader = any_reflection_file('integrated.mtz')
  mtz_object = reader.file_content()
  assert mtz_object.column_labels() == [
    'H', 'K', 'L', 'M_ISYM', 'BATCH', 'IPR', 'SIGIPR', 'I', 'SIGI',
    'FRACTIONCALC', 'XDET', 'YDET', 'ROT', 'LP', 'DQE']
  assert len(mtz_object.batches()) == 45
  batch = mtz_object.batches()[0]
  expected_unit_cell = uctbx.unit_cell((78.07, 78.07, 78.07, 90, 90, 90))
  assert expected_unit_cell.is_similar_to(uctbx.unit_cell(list(batch.cell())))
  assert mtz_object.space_group().type().hall_symbol() == hall_symbol
  # number of reflections we get seems to be rather sensitive to machine
  # or initial conditions?
  assert approx_equal(mtz_object.n_reflections(), 48063, 1e3)
  os.chdir(cwd)


def run(args):
  if not have_dials_regression:
    print "Skipping tst_dials_process.py: dials_regression not available"
    return

  exercises = (exercise_1, exercise_2)
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
