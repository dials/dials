from __future__ import absolute_import, division, print_function

from dials.algorithms.profile_model.gaussian_rs.transform import beam_vector_map

import pytest

@pytest.fixture
def sweep_and_model(dials_regression):
  from dials.model.serialize import load
  import os

  filename = os.path.join(dials_regression, 'centroid_test_data', 'sweep.json')

  class Test(object):
    pass

  storage_class = Test()

  # Load the sweep
  storage_class.sweep = load.sweep(filename)

  # Get the models
  storage_class.beam = storage_class.sweep.get_beam()
  storage_class.detector = storage_class.sweep.get_detector()
  storage_class.gonio = storage_class.sweep.get_goniometer()
  storage_class.scan = storage_class.sweep.get_scan()

  return storage_class

def test_at_corners(sweep_and_model):
  from scitbx import matrix
  from random import randint
  assert len(sweep_and_model.detector) == 1

  # The detector beam vectors
  ds1 = beam_vector_map(sweep_and_model.detector[0], sweep_and_model.beam, True)
  expected_size = sweep_and_model.detector[0].get_image_size()[::-1]
  expected_size = tuple([e + 1 for e in expected_size])
  assert ds1.all() == expected_size

  s0 = sweep_and_model.beam.get_s0()
  m2 = sweep_and_model.gonio.get_rotation_axis()
  s0_length = matrix.col(sweep_and_model.beam.get_s0()).length()

  # Ensure a few random points are correct
  eps = 1e-7
  for k in range(1000):
    j = randint(0, ds1.all()[0]-1)
    i = randint(0, ds1.all()[1]-1)
    y = float(j)
    x = float(i)
    xyz = sweep_and_model.detector[0].get_pixel_lab_coord((x, y))
    s11 = matrix.col(xyz).normalize() * s0_length
    s12 = matrix.col(ds1[j,i])
    assert (s11 - s12).length() <= eps

def test_sub_division_at_corners(sweep_and_model):
  from scitbx import matrix
  from random import randint

  sweep_and_model.n_div = 2

  # The detector beam vectors
  ds1 = beam_vector_map(sweep_and_model.detector[0], sweep_and_model.beam, sweep_and_model.n_div, True)
  expected_size = sweep_and_model.detector[0].get_image_size()[::-1]
  expected_size = tuple([e * sweep_and_model.n_div + 1 for e in expected_size])
  assert ds1.all() == expected_size

  s0 = sweep_and_model.beam.get_s0()
  m2 = sweep_and_model.gonio.get_rotation_axis()
  s0_length = matrix.col(sweep_and_model.beam.get_s0()).length()

  # Ensure a few random points are correct
  eps = 1e-7
  for k in range(1000):
    j = randint(0, ds1.all()[0]-1)
    i = randint(0, ds1.all()[1]-1)
    y = float(j) / sweep_and_model.n_div
    x = float(i) / sweep_and_model.n_div
    xyz = sweep_and_model.detector[0].get_pixel_lab_coord((x, y))
    s11 = matrix.col(xyz).normalize() * s0_length
    s12 = matrix.col(ds1[j,i])
    assert (s11 - s12).length() <= eps

def test_sub_division_at_centres(sweep_and_model):
  from scitbx import matrix
  from random import randint

  sweep_and_model.n_div = 2

  # The detector beam vectors
  ds1 = beam_vector_map(sweep_and_model.detector[0], sweep_and_model.beam, sweep_and_model.n_div, False)
  expected_size = sweep_and_model.detector[0].get_image_size()[::-1]
  expected_size = tuple([e * sweep_and_model.n_div for e in expected_size])
  assert(ds1.all() == expected_size)

  s0 = sweep_and_model.beam.get_s0()
  m2 = sweep_and_model.gonio.get_rotation_axis()
  s0_length = matrix.col(sweep_and_model.beam.get_s0()).length()

  # Ensure a few random points are correct
  eps = 1e-7
  for k in range(1000):
    j = randint(0, ds1.all()[0]-1)
    i = randint(0, ds1.all()[1]-1)
    y = float(j + 0.5) / sweep_and_model.n_div
    x = float(i + 0.5) / sweep_and_model.n_div
    xyz = sweep_and_model.detector[0].get_pixel_lab_coord((x, y))
    s11 = matrix.col(xyz).normalize() * s0_length
    s12 = matrix.col(ds1[j,i])
    assert((s11 - s12).length() <= eps)
