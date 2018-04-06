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
