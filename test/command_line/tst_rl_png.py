from __future__ import absolute_import, division
import os
import libtbx.load_env # required for libtbx.env.find_in_repositories
from libtbx import easy_run
from libtbx.test_utils import open_tmp_directory, approx_equal

def test1():

  dials_regression = libtbx.env.find_in_repositories(
    relative_path="dials_regression",
    test=os.path.isdir)

  data_dir = os.path.join(dials_regression, "centroid_test_data")
  datablock_path = os.path.join(data_dir, "datablock.json")
  strong_pickle = os.path.join(data_dir, "strong.pickle")

  # work in a temporary directory
  cwd = os.path.abspath(os.curdir)
  tmp_dir = open_tmp_directory(suffix="rl_png_test1")
  os.chdir(tmp_dir)

  cmd = 'dials.rl_png %s %s' %(datablock_path, strong_pickle)
  result = easy_run.fully_buffered(command=cmd).raise_if_errors()


  for s in ('beam_vector', 'e3', 'rotation_axis',
            'solution_1', 'solution_2','solution_3'):
    assert os.path.exists('rl_%s.png' %s), s

  os.chdir(cwd)

  return

def test2():

  dials_regression = libtbx.env.find_in_repositories(
    relative_path="dials_regression",
    test=os.path.isdir)

  data_dir = os.path.join(dials_regression, "refinement_test_data", "i04_weak_data")
  datablock_path = os.path.join(data_dir, "experiments.json")
  indexed_pickle = os.path.join(data_dir, "indexed_strong.pickle")

  # work in a temporary directory
  cwd = os.path.abspath(os.curdir)
  tmp_dir = open_tmp_directory(suffix="rl_png_test2")
  os.chdir(tmp_dir)

  cmd = 'dials.rl_png %s %s' %(datablock_path, indexed_pickle)
  result = easy_run.fully_buffered(command=cmd).raise_if_errors()

  for s in ('beam_vector', 'e3', 'rotation_axis', 'a', 'b','c'):
    assert os.path.exists('rl_%s.png' %s), s

  os.chdir(cwd)

  return

def run():
  if not libtbx.env.has_module("dials_regression"):
    print "Skipping tests in " + __file__ + " as dials_regression not present"
    return

  test1()
  test2()

if __name__ == '__main__':
  from dials.test import cd_auto
  with cd_auto(__file__):
    from libtbx.utils import show_times_at_exit
    show_times_at_exit()
    run()
