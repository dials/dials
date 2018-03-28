from __future__ import absolute_import, division, print_function

import os
import math
import random

def test_run(dials_regression, tmpdir):
  tmpdir.chdir()

  from dials.algorithms.profile_model.gaussian_rs import transform
  from dials.algorithms.profile_model.gaussian_rs import CoordinateSystem
  from scitbx import matrix
  from dials.model.serialize import load

  # Set the sweep filename and load the sweep
  filename = os.path.join(dials_regression, 'centroid_test_data', 'sweep.json')
  sweep = load.sweep(filename)

  # Get the models
  beam = sweep.get_beam()
  detector = sweep.get_detector()
  gonio = sweep.get_goniometer()
  assert(len(detector) == 1)

  # Get some stuff
  s0 = beam.get_s0()
  m2 = gonio.get_rotation_axis()
  image_size = detector[0].get_image_size()

  # Get a random s1/phi
  i = random.uniform(0, image_size[0])
  j = random.uniform(1, image_size[1])
  s1 = matrix.col(detector[0].get_pixel_lab_coord((i, j)))
  s1 = s1.normalize() * matrix.col(s0).length()
  phi = random.uniform(0, 5)
  x0 = int(math.floor(i - 10))
  y0 = int(math.floor(j - 10))

  # Set some parameters
  sigma_divergence = beam.get_sigma_divergence(deg=False)
  delta_divergence = 3 * sigma_divergence
  grid_half_size = 4
  step_size = (delta_divergence / grid_half_size,
                    delta_divergence / grid_half_size)

  # Create the coordinate system
  cs = CoordinateSystem(m2, s0, s1, phi)

  # Create the map of s1 coordinates
  s1_map = transform.beam_vector_map(detector[0], beam, True)

  # Create the grid index generator
  generate_indices = transform.GridIndexGenerator(
      cs,
      x0,
      y0,
      step_size,
      grid_half_size,
      s1_map)

  for j in range(0, 20):
    for i in range(0, 20):

      xx = x0 + i
      yy = y0 + j
      if xx < 0 or yy < 0 or xx >= image_size[0] or yy >= image_size[0]:
        continue

      # Get the grid indices
      gi_1, gj_1 = generate_indices(j, i)

      # Get the grid indices
      xyz = matrix.col(detector[0].get_pixel_lab_coord(
          (x0 + i, y0 + j)))
      xyz = xyz.normalize() * matrix.col(s0).length()
      c1, c2 = matrix.col(cs.from_beam_vector(xyz))
      gi_2 = grid_half_size + c1 / step_size[0] + 0.5
      gj_2 = grid_half_size + c2 / step_size[1] + 0.5

      # Check both are the same
      eps = 1e-7
      assert(abs(gj_1 - gj_2) <= eps), (gi_1, gi_2, gj_1, gj_2)
      assert(abs(gi_1 - gi_2) <= eps), (gi_1, gi_2, gj_1, gj_2)
