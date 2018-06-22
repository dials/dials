from __future__ import absolute_import, division, print_function

import os
import pytest
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

  cmd = 'dials.export_bitmaps %s/centroid_0001.cbf prefix=img_ padding=0' %data_dir
  result = easy_run.fully_buffered(cmd).raise_if_errors()
  assert os.path.exists('img_1.png')

  cmd = 'dials.export_bitmaps %s/centroid_0001.cbf prefix=img_ padding=5' %data_dir
  result = easy_run.fully_buffered(cmd).raise_if_errors()
  assert os.path.exists('img_00001.png')

  cmd = 'dials.export_bitmaps %s/centroid_0001.cbf output_file=kittens.png' %data_dir
  result = easy_run.fully_buffered(cmd).raise_if_errors()
  assert os.path.exists('kittens.png')

  with pytest.raises(RuntimeError):
    # setting output filename not allowed with >1 image
    cmd = ' '.join([
      'dials.export_bitmaps', '%s/datablock.json' %data_dir, 'output_file=kittens.png'])
    result = easy_run.fully_buffered(cmd).raise_if_errors()

def test_still_image(dials_regression, tmpdir):
  tmpdir.chdir()
  data_dir = os.path.join(dials_regression, 'image_examples/DLS_I24_stills')

  cmd = 'dials.export_bitmaps %s/still_0001.cbf' %data_dir
  result = easy_run.fully_buffered(cmd).raise_if_errors()

  assert os.path.exists('image0001.png')
