from __future__ import absolute_import, division, print_function

import os
import pytest
from libtbx import easy_run
import glob

def test_export_single_bitmap(dials_regression, run_in_tmpdir):
  data_dir = os.path.join(dials_regression, 'centroid_test_data')

  cmd = 'dials.export_bitmaps %s/centroid_0001.cbf' %data_dir
  result = easy_run.fully_buffered(cmd).raise_if_errors()

  assert os.path.exists('image0001.png')

def test_export_multiple_bitmaps(dials_regression, run_in_tmpdir):
  data_dir = os.path.join(dials_regression, 'centroid_test_data')
  cmd = ' '.join([
    'dials.export_bitmaps', '%s/experiments.json' %data_dir, 'prefix=variance_',
    'binning=2', 'display=variance', 'colour_scheme=inverse_greyscale',
    'brightness=25', 'kernel_size=5,5'])
  result = easy_run.fully_buffered(cmd).raise_if_errors()

  for i in range(1, 8):
    assert os.path.exists('variance_000%i.png' %i)

def test_export_bitmap_with_prefix_and_no_padding(dials_regression, run_in_tmpdir):
  data_dir = os.path.join(dials_regression, 'centroid_test_data')
  cmd = 'dials.export_bitmaps %s/centroid_0001.cbf prefix=img_ padding=0' %data_dir
  result = easy_run.fully_buffered(cmd).raise_if_errors()
  assert os.path.exists('img_1.png')

def test_export_bitmap_with_prefix_and_extra_padding(dials_regression, run_in_tmpdir):
  data_dir = os.path.join(dials_regression, 'centroid_test_data')
  cmd = 'dials.export_bitmaps %s/centroid_0001.cbf prefix=img_ padding=5' %data_dir
  result = easy_run.fully_buffered(cmd).raise_if_errors()
  assert os.path.exists('img_00001.png')

def test_export_bitmap_with_specified_output_filename(dials_regression, run_in_tmpdir):
  data_dir = os.path.join(dials_regression, 'centroid_test_data')
  cmd = 'dials.export_bitmaps %s/centroid_0001.cbf output_file=kittens.png' %data_dir
  result = easy_run.fully_buffered(cmd).raise_if_errors()
  assert os.path.exists('kittens.png')

def test_export_multiple_bitmaps_with_specified_output_filename_fails(dials_regression, run_in_tmpdir):
  data_dir = os.path.join(dials_regression, 'centroid_test_data')
  with pytest.raises(RuntimeError):
    # setting output filename not allowed with >1 image
    cmd = ' '.join([
      'dials.export_bitmaps', '%s/experiments.json' %data_dir, 'output_file=kittens.png'])
    result = easy_run.fully_buffered(cmd).raise_if_errors()

def test_export_still_image(dials_regression, run_in_tmpdir):
  image = os.path.join(dials_regression, 'image_examples', 'DLS_I24_stills', 'still_0001.cbf')

  cmd = 'dials.export_bitmaps %s' % image
  result = easy_run.fully_buffered(cmd).raise_if_errors()

  assert os.path.exists('image0001.png')

def test_export_multi_panel(dials_regression, run_in_tmpdir):
  image = os.path.join(dials_regression, 'image_examples', 'DLS_I23', 'germ_13KeV_0001.cbf')

  for binning in (1, 4):
    cmd = 'dials.export_bitmaps %s binning=%i prefix=binning_%i_' % (
      image, binning, binning)
    result = easy_run.fully_buffered(cmd).raise_if_errors()

    assert os.path.exists('binning_%i_0001.png' % binning)

def test_export_restricted_multiimage(dials_regression, run_in_tmpdir):
  "Test exporting a subset of an imageset"
  data_dir = os.path.join(dials_regression, 'centroid_test_data')

  cmd = ' '.join([
    'dials.export_bitmaps', '%s/experiments.json' %data_dir, 'imageset_index=2'])
  easy_run.fully_buffered(cmd).raise_if_errors()
  assert glob.glob("*.png") == ["image0002.png"], "Only one image exported"
