from __future__ import absolute_import, division, print_function

import cPickle as pickle
from glob import glob
import os

import procrunner
import pytest

from dials.array_family import flex # import dependency

def test_find_spots_from_images(dials_regression, tmpdir):
  tmpdir.chdir()

  result = procrunner.run_process([
      "dials.find_spots",
      "output.reflections=spotfinder.pickle",
      "output.shoeboxes=True",
  ] + glob(os.path.join(dials_regression, "centroid_test_data", "centroid*.cbf")))
  assert result['exitcode'] == 0
  assert result['stderr'] == ''
  assert os.path.exists("spotfinder.pickle")

  with open("spotfinder.pickle", "rb") as f:
    reflections = pickle.load(f)
  assert len(reflections) == 653
  refl = reflections[0]
  assert refl['intensity.sum.value'] == pytest.approx(42)
  assert refl['bbox'] == pytest.approx((1398, 1400, 513, 515, 0, 1))
  assert refl['xyzobs.px.value'] == pytest.approx((1399.1190476190477, 514.2142857142857, 0.5))
  assert "shoebox" in reflections

def test_find_spots_with_resolution_filter(dials_regression, tmpdir):
  tmpdir.chdir()

  result = procrunner.run_process([
      "dials.find_spots",
      "output.reflections=spotfinder.pickle",
      "output.shoeboxes=False",
      "filter.d_min=2",
      "filter.d_max=15",
  ] + glob(os.path.join(dials_regression, "centroid_test_data", "centroid*.cbf")))
  assert result['exitcode'] == 0
  assert result['stderr'] == ''
  assert os.path.exists("spotfinder.pickle")

  with open("spotfinder.pickle", "rb") as f:
    reflections = pickle.load(f)
  assert len(reflections) == 467
  assert "shoebox" not in reflections

def test_find_spots_with_hot_mask(dials_regression, tmpdir):
  tmpdir.chdir()

  # now write a hot mask
  result = procrunner.run_process([
      "dials.find_spots",
      "write_hot_mask=True",
      "output.reflections=spotfinder.pickle",
      "output.shoeboxes=False",
  ] + glob(os.path.join(dials_regression, "centroid_test_data", "centroid*.cbf")))
  assert result['exitcode'] == 0
  assert result['stderr'] == ''
  assert os.path.exists("spotfinder.pickle")
  assert os.path.exists("hot_mask_0.pickle")

  with open("spotfinder.pickle", "rb") as f:
    reflections = pickle.load(f)
  assert len(reflections) == 653
  assert "shoebox" not in reflections

  with open("hot_mask_0.pickle", "rb") as f:
    mask = pickle.load(f)
  assert len(mask) == 1
  assert mask[0].count(False) == 12

def test_find_spots_with_hot_mask_with_prefix(dials_regression, tmpdir):
  tmpdir.chdir()

  # now write a hot mask
  result = procrunner.run_process([
      "dials.find_spots",
      "write_hot_mask=True",
      "hot_mask_prefix=my_hot_mask",
      "output.reflections=spotfinder.pickle",
      "output.shoeboxes=False",
  ] + glob(os.path.join(dials_regression, "centroid_test_data", "centroid*.cbf")))
  assert result['exitcode'] == 0
  assert result['stderr'] == ''
  assert os.path.exists("spotfinder.pickle")
  assert os.path.exists("my_hot_mask_0.pickle")

  with open("spotfinder.pickle", "rb") as f:
    reflections = pickle.load(f)
  assert len(reflections) == 653
  assert "shoebox" not in reflections
  with open("my_hot_mask_0.pickle", "rb") as f:
    mask = pickle.load(f)
  assert len(mask) == 1
  assert mask[0].count(False) == 12

def test_find_spots_with_generous_parameters(dials_regression, tmpdir):
  tmpdir.chdir()

  # now with more generous parameters
  result = procrunner.run_process([
      "dials.find_spots",
      "min_spot_size=3",
      "max_separation=3",
      "output.reflections=spotfinder.pickle",
  ] + glob(os.path.join(dials_regression, "centroid_test_data", "centroid*.cbf")))
  assert result['exitcode'] == 0
  assert result['stderr'] == ''
  assert os.path.exists("spotfinder.pickle")

  with open("spotfinder.pickle", "rb") as f:
    reflections = pickle.load(f)
  assert len(reflections) == 678

def test_find_spots_with_user_defined_mask(dials_regression, tmpdir):
  tmpdir.chdir()

  # Now with a user defined mask
  result = procrunner.run_process([
      "dials.find_spots",
      "output.reflections=spotfinder.pickle",
      "output.shoeboxes=True",
      "lookup.mask=" + os.path.join(dials_regression, "centroid_test_data", "mask.pickle"),
  ] + glob(os.path.join(dials_regression, "centroid_test_data", "centroid*.cbf")))
  assert result['exitcode'] == 0
  assert result['stderr'] == ''
  assert os.path.exists("spotfinder.pickle")

  with open("spotfinder.pickle", "rb") as f:
    reflections = pickle.load(f)

  from dxtbx.datablock import DataBlockFactory
  datablocks = DataBlockFactory.from_json_file(
      os.path.join(dials_regression, "centroid_test_data", "datablock.json"))
  assert len(datablocks) == 1
  imageset = datablocks[0].extract_imagesets()[0]
  detector = imageset.get_detector()
  beam = imageset.get_beam()
  for x, y, z in reflections['xyzobs.px.value']:
    d = detector[0].get_resolution_at_pixel(beam.get_s0(), (x, y))
    assert d >= 3

def test_find_spots_with_user_defined_region(dials_regression, tmpdir):
  tmpdir.chdir()

  result = procrunner.run_process([
      "dials.find_spots",
      "output.reflections=spotfinder.pickle",
      "output.shoeboxes=True",
      "region_of_interest=800,1200,800,1200",
  ] + glob(os.path.join(dials_regression, "centroid_test_data", "centroid*.cbf")))
  assert result['exitcode'] == 0
  assert result['stderr'] == ''
  assert os.path.exists("spotfinder.pickle")

  with open("spotfinder.pickle", "rb") as f:
    reflections = pickle.load(f)
  x, y, z = reflections['xyzobs.px.value'].parts()
  assert x.all_ge(800)
  assert y.all_ge(800)
  assert x.all_lt(1200)
  assert y.all_lt(1200)

def test_find_spots_with_xfel_stills(dials_regression, tmpdir):
  tmpdir.chdir()

  # now with XFEL stills
  result = procrunner.run_process([
      "dials.find_spots",
      os.path.join(dials_regression, "spotfinding_test_data", "idx-s00-20131106040302615.cbf"),
      "output.reflections=spotfinder.pickle",
  ])
  assert result['exitcode'] == 0
  assert result['stderr'] == ''
  assert os.path.exists("spotfinder.pickle")

  with open("spotfinder.pickle", "rb") as f:
    reflections = pickle.load(f)
  assert len(reflections) == 2643

def test_polygon():
  from dials.algorithms.polygon import polygon

  x = 1
  y = 1
  vertices = [(0,0), (2,0), (2,2), (0,2)]
  poly = polygon(vertices)
  assert poly.is_inside(x,y)

  poly = polygon([(3,5), (40,90), (80,70), (50,50), (70,20)])
  for p in [(42,40), (16,25), (64,67), (16,30), (48, 45), (21,30)]:
    assert poly.is_inside(p[0], p[1])
  for p in [(59,4), (70,15), (57,14), (21,78), (37,100), (88,89)]:
    assert not poly.is_inside(p[0], p[1])

  if 0:
    # for visual confirmation of algorithm
    from scitbx.array_family import flex
    inside_points = flex.vec2_double()
    outside_points = flex.vec2_double()
    import random
    x_max = 100
    y_max = 100
    for i in range(1000):
      x = random.randint(0,x_max)
      y = random.randint(0,y_max)
      is_inside = poly.is_inside(x, y)
      if is_inside:
        inside_points.append((x,y))
      else:
        outside_points.append((x,y))
    from matplotlib import pyplot
    from matplotlib.patches import Polygon
    v = poly.vertices + poly.vertices[:1]
    fig = pyplot.figure()
    ax = fig.add_subplot(111)
    ax.add_patch(Polygon(poly.vertices, closed=True, fill=False))
    inside_x, inside_y = inside_points.parts()
    outside_x, outside_y = outside_points.parts()
    ax.scatter(inside_x, inside_y, marker='+', c='r')
    ax.scatter(outside_x, outside_y, marker='+', c='b')
    pyplot.show()
