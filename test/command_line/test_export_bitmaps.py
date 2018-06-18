from __future__ import absolute_import, division, print_function

import os
from libtbx import easy_run

def test_export_bitmaps(dials_regression, tmpdir):
  tmpdir.chdir()
  data_dir = os.path.join(dials_regression, 'centroid_test_data')

  cmd = 'dials.export_bitmaps %s/centroid_0001.cbf' %data_dir
  result = easy_run.fully_buffered(cmd).raise_if_errors()

  assert os.path.exists('image0001.png')

  cmd = ' '.join([
    'dials.export_bitmaps', '%s/datablock.json' %data_dir, 'prefix=variance_',
    'binning=2', 'display=variance', 'colour_scheme=inverse_greyscale',
    'brightness=25', 'kernel_size=5,5'])
  result = easy_run.fully_buffered(cmd).raise_if_errors()

  for i in range(1, 8):
    assert os.path.exists('variance_000%i.png' %i)

def test_still_image(dials_regression, tmpdir):
  tmpdir.chdir()
  data_dir = os.path.join(dials_regression, 'image_examples/DLS_I24_stills')

  cmd = 'dials.export_bitmaps %s/still_0001.cbf' %data_dir
  result = easy_run.fully_buffered(cmd).raise_if_errors()

  assert os.path.exists('image0001.png')
