from __future__ import absolute_import, division, print_function

import os

def test_polygon():
  from dials.util.ext import is_inside_polygon
  from scitbx.array_family import flex

  poly = flex.vec2_double(((0,0), (1,0), (1,1), (0,1)))
  assert is_inside_polygon(poly, 0, 0)
  assert is_inside_polygon(poly, 0.5, 0.5)
  assert not is_inside_polygon(poly, 1, 1.01)

  points = flex.vec2_double(((0.3, 0.8), (0.3, 1.5), (-8,9), (0.00001, 0.9999)))
  assert list(is_inside_polygon(poly, points)) == [True, False, False, True]

def test_dynamic_shadowing(dials_regression):
  import libtbx

  def exercise_one_image(
    path, count_only_shadow, count_mask_shadow, count_mask_no_shadow):

    from dxtbx.datablock import DataBlockFactory
    assert os.path.exists(path), path
    for shadowing in (libtbx.Auto, True, False):
      format_kwargs = {'dynamic_shadowing': shadowing}
      datablock = DataBlockFactory.from_filenames([path], format_kwargs=format_kwargs)[0]
      imageset = datablock.extract_imagesets()[0]
      detector = imageset.get_detector()
      scan = imageset.get_scan()
      filename = imageset.get_path(0)
      masker = imageset.masker().format_class(filename, **format_kwargs).get_goniometer_shadow_masker()
      assert masker is not None
      mask = masker.get_mask(detector, scan_angle=scan.get_oscillation()[0])
      assert len(mask) == len(detector)
      # only shadowed pixels masked
      assert mask[0].count(False) == count_only_shadow, (mask[0].count(False), count_only_shadow)
      mask = imageset.get_mask(0)
      # dead pixels, pixels in gaps, etc also masked
      if shadowing is libtbx.Auto or shadowing is True:
        assert mask[0].count(False) == count_mask_shadow, (mask[0].count(False), count_mask_shadow)
      else:
        assert mask[0].count(False) == count_mask_no_shadow, (mask[0].count(False), count_mask_no_shadow)

  exercise_one_image(
    path=os.path.join(dials_regression,
      "shadow_test_data/DLS_I04_SmarGon/Th_3_O45_C45_P48_1_0500.cbf"),
    count_only_shadow=426758,
    count_mask_shadow=917940,
    count_mask_no_shadow=528032
  )
  exercise_one_image(
    path=os.path.join(dials_regression,
      "shadow_test_data/DLS_I03_SmarGon/protk_2_0600.cbf"),
    count_only_shadow=519100,
    count_mask_shadow=1002068,
    count_mask_no_shadow=527314
  )

def test_shadow_plot(dials_regression, run_in_tmpdir):
  path = os.path.join(
    dials_regression, "shadow_test_data/DLS_I04_SmarGon/Th_3_O45_C45_P48_1_0500.cbf")

  from libtbx.easy_run import fully_buffered
  from libtbx.test_utils import approx_equal

  result = fully_buffered('dials.import %s' %path).raise_if_errors()
  result = fully_buffered(
    'dials.shadow_plot datablock.json json=shadow.json').raise_if_errors()
  assert os.path.exists('scan_shadow_plot.png')
  assert os.path.exists('shadow.json')
  with open('shadow.json', 'rb') as f:
    import json
    d = json.load(f)
    assert d.keys() == ['fraction_shadowed', 'scan_points']
    assert approx_equal(d['fraction_shadowed'], [0.06856597327776767])
    assert approx_equal(d['scan_points'], [94.9])

  result = fully_buffered(
    'dials.shadow_plot datablock.json mode=2d plot=shadow_2d.png').raise_if_errors()
  assert os.path.exists('shadow_2d.png')

def test_filter_shadowed_reflections(dials_regression):
  experiments_json = os.path.join(
    dials_regression, "shadow_test_data/DLS_I04_SmarGon/experiments.json")

  predicted_pickle = os.path.join(
    dials_regression, "shadow_test_data/DLS_I04_SmarGon/predicted.pickle")

  from dxtbx.serialize import load
  experiments = load.experiment_list(experiments_json, check_format=True)

  from libtbx import easy_pickle
  predicted = easy_pickle.load(predicted_pickle)
  from dials.algorithms.shadowing.filter import filter_shadowed_reflections
  for experiment_goniometer in (True, False):
    shadowed = filter_shadowed_reflections(
      experiments, predicted, experiment_goniometer=experiment_goniometer)
    assert shadowed.count(True) == 17
    assert shadowed.count(False) == 674
