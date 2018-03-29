from __future__ import absolute_import, division

import os

import pytest

@pytest.fixture
def setup(dials_regression):
  from dials.algorithms.profile_model.gaussian_rs import BBoxCalculator3D
  from dials.model.serialize import load
  from math import pi

  filename = os.path.join(dials_regression,
      'centroid_test_data', 'sweep.json')

  sweep = load.sweep(filename)

  fixture = {}

  # Get the models
  fixture['beam'] = sweep.get_beam()
  fixture['detector'] = sweep.get_detector()
  fixture['gonio'] = sweep.get_goniometer()
  fixture['scan'] = sweep.get_scan()

  # Set the delta_divergence/mosaicity
  n_sigma = 5
  sigma_divergence = 0.060 * pi / 180
  mosaicity = 0.154 * pi / 180
  fixture['delta_divergence'] = n_sigma * sigma_divergence
  fixture['delta_mosaicity'] = n_sigma * mosaicity

  # Create the bounding box calculator
  fixture['calculate_bbox'] = BBoxCalculator3D(
      fixture['beam'], fixture['detector'],
      fixture['gonio'], fixture['scan'],
      fixture['delta_divergence'],
      fixture['delta_mosaicity'])
  return fixture

def test_outer_bounds(setup):
    from scitbx import matrix
    from random import uniform
    from dials.algorithms.profile_model.gaussian_rs import CoordinateSystem

    assert(len(setup['detector']) == 1)
    s0 = setup['beam'].get_s0()
    m2 = setup['gonio'].get_rotation_axis()
    s0_length = matrix.col(setup['beam'].get_s0()).length()
    width, height = setup['detector'][0].get_image_size()
    zrange = setup['scan'].get_array_range()

    for i in range(1000):

      # Get random x, y, z
      x = uniform(0, 2000)
      y = uniform(0, 2000)
      z = uniform(0, 9)

      # Get random s1, phi, panel
      s1 = matrix.col(setup['detector'][0].get_pixel_lab_coord(
          (x, y))).normalize() * s0_length
      phi = setup['scan'].get_angle_from_array_index(z, deg=False)
      panel = 0

      # Calculate the bounding box
      bbox = setup['calculate_bbox'](s1, z, panel)
      x1, x2 = bbox[0], bbox[1]
      y1, y2 = bbox[2], bbox[3]
      z1, z2 = bbox[4], bbox[5]

      # Calculate the rotation angle for each point
      phi_dash1 = setup['scan'].get_angle_from_array_index(z1, deg = False)
      phi_dash2 = setup['scan'].get_angle_from_array_index(z2, deg = False)

      # Create the XDS coordinate system
      xcs = CoordinateSystem(m2, s0, s1, phi)

      # Calculate reciprocal space coordinates at each point
      e11, e21, e31 = xcs.from_beam_vector_and_rotation_angle(s1, phi_dash1)
      e12, e22, e32 = xcs.from_beam_vector_and_rotation_angle(s1, phi_dash2)

      # Check vertical edges
      for j in range(bbox[2], bbox[3] + 1):
        xyz1 = setup['detector'][0].get_pixel_lab_coord((bbox[0], j))
        xyz2 = setup['detector'][0].get_pixel_lab_coord((bbox[1] + 1, j))
        sdash1 = matrix.col(xyz1).normalize() * s0_length
        sdash2 = matrix.col(xyz2).normalize() * s0_length
        e11, e21, e3 = xcs.from_beam_vector_and_rotation_angle(sdash1, phi)
        e12, e22, e3 = xcs.from_beam_vector_and_rotation_angle(sdash2, phi)
        if bbox[0] > 0 and bbox[1] < width:
          assert(abs(e11) >= setup['delta_divergence'] or
                 abs(e21) >= setup['delta_divergence'])
          assert(abs(e12) >= setup['delta_divergence'] or
                 abs(e22) >= setup['delta_divergence'])

      # Check horizontal edges
      for i in range(bbox[0], bbox[1] + 1):
        xyz1 = setup['detector'][0].get_pixel_lab_coord((i, bbox[2]))
        xyz2 = setup['detector'][0].get_pixel_lab_coord((i, bbox[3] + 1))
        sdash1 = matrix.col(xyz1).normalize() * s0_length
        sdash2 = matrix.col(xyz2).normalize() * s0_length
        e11, e21, e3 = xcs.from_beam_vector_and_rotation_angle(sdash1, phi)
        e12, e22, e3 = xcs.from_beam_vector_and_rotation_angle(sdash2, phi)
        if bbox[2] > 0 and bbox[3] < height:
          assert(abs(e11) >= setup['delta_divergence'] or
                 abs(e21) >= setup['delta_divergence'])
          assert(abs(e12) >= setup['delta_divergence'] or
                 abs(e22) >= setup['delta_divergence'])

      # All e3 coords >= delta_mosaicity
      if bbox[4] > zrange[0] and bbox[5] < zrange[1]:
        assert(abs(e31) >= setup['delta_mosaicity'])
        assert(abs(e32) >= setup['delta_mosaicity'])


def test_radius(setup):
    from scitbx import matrix
    from random import uniform
    from dials.algorithms.profile_model.gaussian_rs import CoordinateSystem
    from math import sqrt

    s0 = setup['beam'].get_s0()
    m2 = setup['gonio'].get_rotation_axis()
    s0_length = matrix.col(setup['beam'].get_s0()).length()

    width, height = setup['detector'][0].get_image_size()
    zrange = setup['scan'].get_array_range()
    radius12 = setup['delta_divergence']
    radius3 = setup['delta_mosaicity']

    for i in range(1000):

      # Get random x, y, z
      x = uniform(0, 2000)
      y = uniform(0, 2000)
      z = uniform(0, 9)

      # Get random s1, phi, panel
      s1 = matrix.col(setup['detector'][0].get_pixel_lab_coord(
          (x, y))).normalize() * s0_length
      phi = setup['scan'].get_angle_from_array_index(z, deg=False)
      panel = 0

      # Calculate the bounding box
      bbox = setup['calculate_bbox'](s1, z, panel)
      x1, x2 = bbox[0], bbox[1]
      y1, y2 = bbox[2], bbox[3]
      z1, z2 = bbox[4], bbox[5]

      # Calculate the rotation angle for each point
      phi_dash1 = setup['scan'].get_angle_from_array_index(z1, deg = False)
      phi_dash2 = setup['scan'].get_angle_from_array_index(z2, deg = False)

      # Create the XDS coordinate system
      xcs = CoordinateSystem(m2, s0, s1, phi)

      # Check vertical edges
      for j in range(bbox[2], bbox[3] + 1):
        xyz1 = setup['detector'][0].get_pixel_lab_coord((bbox[0], j))
        xyz2 = setup['detector'][0].get_pixel_lab_coord((bbox[1] + 1, j))
        sdash1 = matrix.col(xyz1).normalize() * s0_length
        sdash2 = matrix.col(xyz2).normalize() * s0_length
        e11, e21, e31 = xcs.from_beam_vector_and_rotation_angle(sdash1, phi)
        e12, e22, e31 = xcs.from_beam_vector_and_rotation_angle(sdash2, phi)
        if bbox[0] > 0 and bbox[1] < width:
          assert(sqrt(e11**2 + e21**2) >= radius12)
          assert(sqrt(e12**2 + e22**2) >= radius12)

      # Check horizontal edges
      for i in range(bbox[0], bbox[1] + 1):
        xyz1 = setup['detector'][0].get_pixel_lab_coord((i, bbox[2]))
        xyz2 = setup['detector'][0].get_pixel_lab_coord((i, bbox[3] + 1))
        sdash1 = matrix.col(xyz1).normalize() * s0_length
        sdash2 = matrix.col(xyz2).normalize() * s0_length
        e11, e21, e32 = xcs.from_beam_vector_and_rotation_angle(sdash1, phi)
        e12, e22, e32 = xcs.from_beam_vector_and_rotation_angle(sdash2, phi)
        if bbox[2] > 0 and bbox[3] < height:
          assert(sqrt(e11**2 + e21**2) >= radius12)
          assert(sqrt(e12**2 + e22**2) >= radius12)
