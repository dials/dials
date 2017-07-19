from __future__ import absolute_import, division

def exercise_polygon():
  from dials.util import is_inside_polygon
  from scitbx.array_family import flex

  poly = flex.vec2_double(((0,0), (1,0), (1,1), (0,1)))
  assert is_inside_polygon(poly, 0, 0)
  assert is_inside_polygon(poly, 0.5, 0.5)
  assert not is_inside_polygon(poly, 1, 1.01)

  points = flex.vec2_double(((0.3, 0.8), (0.3, 1.5), (-8,9), (0.00001, 0.9999)))
  assert list(is_inside_polygon(poly, points)) == [True, False, False, True]

def exercise_dynamic_shadowing():
  import os
  import libtbx.load_env
  try:
    dials_regression = libtbx.env.dist_path('dials_regression')
  except KeyError, e:
    print 'SKIP: dials_regression not configured'
    exit(0)

  path = os.path.join(
    dials_regression, "shadow_test_data/DLS_I04_SmarGon/Th_3_O45_C45_P48_1_0500.cbf")

  from dxtbx.datablock import DataBlockFactory
  format_kwargs = {'dynamic_shadowing': True}
  assert os.path.exists(path), path
  datablock = DataBlockFactory.from_filenames([path], format_kwargs=format_kwargs)[0]
  imageset = datablock.extract_imagesets()[0]
  detector = imageset.get_detector()
  scan = imageset.get_scan()
  masker = imageset.reader().get_format().get_goniometer_shadow_masker()
  assert masker is not None
  mask = masker.get_mask(detector, scan_angle=scan.get_oscillation()[0])
  assert len(mask) == len(detector)
  # only shadowed pixels masked
  assert (mask[0].count(True), mask[0].count(False)) == (5797243, 426758)
  mask = imageset.get_mask(0)
  # dead pixels, pixels in gaps, etc also masked
  assert (mask[0].count(True), mask[0].count(False)) == (5306061, 917940)


def run():
  exercise_polygon()
  exercise_dynamic_shadowing()
  print 'OK'

if __name__ == '__main__':
  run()
