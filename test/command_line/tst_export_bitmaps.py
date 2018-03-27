from __future__ import absolute_import, division

import os

import libtbx.load_env
from libtbx import easy_run

dials_regression = libtbx.env.find_in_repositories('dials_regression')

def exercise_export_bitmaps():
  if dials_regression is None:
    print 'Skipping exercise_export_bitmaps(): dials_regression not available'
    return

  data_dir = os.path.join(dials_regression, 'centroid_test_data')

  cmd = 'dials.export_bitmaps %s/centroid_0001.cbf' %data_dir
  print cmd
  result = easy_run.fully_buffered(cmd).raise_if_errors()

  assert os.path.exists('image0001.png')

  cmd = ' '.join([
    'dials.export_bitmaps', '%s/datablock.json' %data_dir, 'prefix=variance_',
    'binning=2', 'display=variance', 'colour_scheme=inverse_greyscale',
    'brightness=25', 'kernel_size=5,5'])
  print cmd
  result = easy_run.fully_buffered(cmd).raise_if_errors()

  for i in range(1, 8):
    assert os.path.exists('variance_000%i.png' %i)

def run():
  exercise_export_bitmaps()

if __name__ == '__main__':
  from dials.test import cd_auto
  with cd_auto(__file__):
    run()
