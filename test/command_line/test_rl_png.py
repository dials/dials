from __future__ import absolute_import, division, print_function

import os

from libtbx import easy_run

def test_rl_png_datablock(dials_regression, tmpdir):
  tmpdir.chdir()

  data_dir = os.path.join(dials_regression, "centroid_test_data")
  datablock_path = os.path.join(data_dir, "datablock.json")
  strong_pickle = os.path.join(data_dir, "strong.pickle")

  cmd = 'dials.rl_png %s %s' %(datablock_path, strong_pickle)
  result = easy_run.fully_buffered(command=cmd).raise_if_errors()

  for s in ('beam_vector', 'e3', 'rotation_axis',
            'solution_1', 'solution_2','solution_3'):
    assert os.path.exists('rl_%s.png' %s), s


def test_rl_png_experiments(dials_regression, tmpdir):
  tmpdir.chdir()

  data_dir = os.path.join(dials_regression, "refinement_test_data", "i04_weak_data")
  datablock_path = os.path.join(data_dir, "experiments.json")
  indexed_pickle = os.path.join(data_dir, "indexed_strong.pickle")

  cmd = 'dials.rl_png %s %s' %(datablock_path, indexed_pickle)
  result = easy_run.fully_buffered(command=cmd).raise_if_errors()

  for s in ('beam_vector', 'e3', 'rotation_axis', 'a', 'b','c'):
    assert os.path.exists('rl_%s.png' %s), s
