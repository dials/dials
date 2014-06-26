from __future__ import division
import glob
import os
from libtbx import easy_run
from libtbx.test_utils import approx_equal, open_tmp_directory
from cctbx import uctbx

def exercise_test_case_1():
  data_dir = os.environ.get('DIALS_PROCESS_TEST_CASE_1_DATA_DIR')
  if data_dir is None:
    print "Skipping exercise_test_case_1(): DIALS_PROCESS_TEST_CASE_1_DATA_DIR not set"
    return
  cwd = os.path.abspath(os.curdir)
  tmp_dir = open_tmp_directory()
  os.chdir(tmp_dir)
  print tmp_dir
  g = glob.glob(os.path.join(data_dir, "PNAS_M3S1_2_*cbf"))
  cmd = " ".join(["dials.process", "scan_varying=True"]
                + ['"%s"' %p for p in g]
                 )
  #print cmd
  result = easy_run.fully_buffered(cmd).raise_if_errors()
  #result = easy_run.call(cmd)
  for out_file in ['datablock.json', 'experiments.json', 'extracted.tar',
                   'indexed.pickle', 'integrated.mtz', 'integrated.pickle',
                   'peaks.pdb', 'strong.pickle']:
    assert os.path.exists(out_file)

  from iotbx.reflection_file_reader import any_reflection_file
  reader = any_reflection_file('integrated.mtz')
  mtz_object = reader.file_content()
  assert mtz_object.n_reflections() == 188338
  assert mtz_object.column_labels() == [
    'H', 'K', 'L', 'M_ISYM', 'BATCH', 'I', 'SIGI', 'FRACTIONCALC',
    'XDET', 'YDET', 'ROT']
  assert len(mtz_object.batches()) == 1800
  batch = mtz_object.batches()[0]
  expected_unit_cell = uctbx.unit_cell(
      [39.385, 42.30, 42.30, 90.0, 90.0, 90.0])
  assert expected_unit_cell.is_similar_to(uctbx.unit_cell(list(batch.cell())))
  os.chdir(cwd)

def run():
  exercise_test_case_1()
  print "OK"

if __name__ == '__main__':
  run()
