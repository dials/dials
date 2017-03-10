from __future__ import absolute_import, division
import os
import libtbx.load_env
from libtbx import easy_run
from libtbx.test_utils import approx_equal
from libtbx.test_utils import open_tmp_directory
from cctbx import uctbx


have_dials_regression = libtbx.env.has_module("dials_regression")
if have_dials_regression:
  dials_regression = libtbx.env.find_in_repositories(
    relative_path="dials_regression",
    test=os.path.isdir)


def exercise_refine_bravais_settings():
  if not have_dials_regression:
    print "Skipping exercise_refine_bravais_settings(): dials_regression not available."
    return

  data_dir = os.path.join(dials_regression, "indexing_test_data", "i04_weak_data")
  pickle_path = os.path.join(data_dir, "indexed.pickle")
  experiments_path = os.path.join(data_dir, "experiments.json")
  commands = ["dials.refine_bravais_settings",
              pickle_path,
              experiments_path,
              "reflections_per_degree=5",
              "minimum_sample_size=500",
              "beam.fix=all",
              "detector.fix=all",
              "prefix=tst_"]
  command = " ".join(commands)
  print command
  cwd = os.path.abspath(os.curdir)
  tmp_dir = open_tmp_directory()
  os.chdir(tmp_dir)
  result = easy_run.fully_buffered(command=command).raise_if_errors()
  for i in range(1, 10):
    assert os.path.exists("tst_bravais_setting_%i.json" %i)
  from dxtbx.serialize import load
  experiments_list = load.experiment_list(
    "tst_bravais_setting_9.json", check_format=False)
  assert len(experiments_list) == 1
  assert experiments_list[0].crystal.get_unit_cell().is_similar_to(
    uctbx.unit_cell((57.782, 57.782, 150.011, 90, 90, 90)))
  assert experiments_list[0].crystal.get_space_group().type().hall_symbol() \
         == ' P 4'

  assert os.path.exists("tst_bravais_summary.json")
  from json import load
  bravais_summary = load(open("tst_bravais_summary.json", "rb"))
  assert bravais_summary.keys() == [
    '1', '3', '2', '5', '4', '7', '6', '9', '8']
  bravais_summary['9'].keys() == [
    'bravais', 'max_angular_difference', 'unit_cell', 'rmsd', 'nspots']

  assert approx_equal(
    bravais_summary['9']['unit_cell'],
    [57.78, 57.78, 150.0, 90.0, 90.0, 90.0], eps=1e-1)
  assert bravais_summary['9']['bravais'] == 'tP'
  assert bravais_summary['9']['recommended'] == True
  assert approx_equal(bravais_summary['9']['rmsd'], 0.047, eps=1e-2)
  os.chdir(cwd)


def exercise_refine_bravais_settings_2():
  if not have_dials_regression:
    print "Skipping exercise_refine_bravais_settings(): dials_regression not available."
    return

  data_dir = os.path.join(dials_regression, "indexing_test_data", "multi_sweep")
  pickle_path = os.path.join(data_dir, "indexed.pickle")
  experiments_path = os.path.join(data_dir, "experiments.json")
  commands = ["dials.refine_bravais_settings", pickle_path, experiments_path]
  command = " ".join(commands)
  print command
  cwd = os.path.abspath(os.curdir)
  tmp_dir = open_tmp_directory()
  os.chdir(tmp_dir)
  result = easy_run.fully_buffered(command=command).raise_if_errors()
  for i in range(1, 10):
    assert os.path.exists("bravais_setting_%i.json" %i)
  from dxtbx.serialize import load
  experiments_list = load.experiment_list(
    "bravais_setting_9.json", check_format=False)
  assert len(experiments_list) == 4
  assert len(experiments_list.crystals()) == 1
  assert experiments_list[0].crystal.get_unit_cell().is_similar_to(
    uctbx.unit_cell((7.31, 7.31, 6.82, 90.00, 90.00, 90.00)))
  assert experiments_list[0].crystal.get_space_group().type().hall_symbol() \
         == ' I 4'
  assert os.path.exists("bravais_summary.json")
  from json import load
  bravais_summary = load(open("bravais_summary.json", "rb"))
  for i in range(1, 23): assert str(i) in bravais_summary.keys()

  assert approx_equal(
    bravais_summary['9']['unit_cell'],
    [7.31, 7.31, 6.82, 90.00, 90.00, 90.00], eps=1e-1)
  assert bravais_summary['9']['bravais'] == 'tI'
  assert approx_equal(bravais_summary['9']['rmsd'], 0.103, eps=1e-2)
  assert bravais_summary['9']['recommended'] == True
  os.chdir(cwd)

def exercise_refine_bravais_settings_3():
  if not have_dials_regression:
    print "Skipping exercise_refine_bravais_settings(): dials_regression not available."
    return

  data_dir = os.path.join(dials_regression, "indexing_test_data", "trypsin")
  pickle_path = os.path.join(data_dir, "indexed.pickle")
  experiments_path = os.path.join(data_dir, "experiments.json")
  commands = ["dials.refine_bravais_settings",
              pickle_path,
              experiments_path,
              "crystal_id=1"]
  command = " ".join(commands)
  print command
  cwd = os.path.abspath(os.curdir)
  tmp_dir = open_tmp_directory()
  os.chdir(tmp_dir)
  result = easy_run.fully_buffered(command=command).raise_if_errors()
  for i in range(1, 10):
    assert os.path.exists("bravais_setting_%i.json" %i)
  from dxtbx.serialize import load
  experiments_list = load.experiment_list(
    "bravais_setting_5.json", check_format=False)
  assert len(experiments_list) == 1
  assert experiments_list[0].crystal.get_unit_cell().is_similar_to(
    uctbx.unit_cell((54.37, 58.29, 66.51, 90.00, 90.00, 90.00)))
  assert experiments_list[0].crystal.get_space_group().type().hall_symbol() \
         == ' P 2 2'

  assert os.path.exists("bravais_summary.json")
  from json import load
  bravais_summary = load(open("bravais_summary.json", "rb"))
  assert bravais_summary.keys() == [
    '1', '3', '2', '5', '4', '7', '6', '9', '8']

  assert approx_equal(
    bravais_summary['5']['unit_cell'],
    [54.37, 58.29, 66.51, 90.00, 90.00, 90.00], eps=1e-1)
  assert bravais_summary['5']['bravais'] == 'oP'
  assert approx_equal(bravais_summary['5']['rmsd'], 0.1200, eps=1e-2)
  assert bravais_summary['5']['recommended'] == True
  assert bravais_summary['9']['recommended'] == False
  os.chdir(cwd)



def run():
  exercise_refine_bravais_settings()
  exercise_refine_bravais_settings_2()
  exercise_refine_bravais_settings_3()

if __name__ == '__main__':
  from dials.test import cd_auto
  with cd_auto(__file__):
    run()
    print "OK"
