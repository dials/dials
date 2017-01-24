from __future__ import absolute_import, division

class Test(object):

  def __init__(self):
    from dials.algorithms.profile_model.gaussian_rs import transform
    from dials.algorithms.profile_model.gaussian_rs import CoordinateSystem
    import os
    import libtbx.load_env
    from math import floor
    from scitbx import matrix
    from dials.model.serialize import load
    from random import uniform

    try:
      dials_regression = libtbx.env.dist_path('dials_regression')
    except KeyError, e:
      print 'FAIL: dials_regression not configured'
      exit(0)

    # Set the sweep filename and load the sweep
    filename = os.path.join(
        dials_regression,
        'centroid_test_data',
        'sweep.json')
    self.sweep = load.sweep(filename)

    # Get the models
    self.beam = self.sweep.get_beam()
    self.detector = self.sweep.get_detector()
    self.gonio = self.sweep.get_goniometer()
    assert(len(self.detector) == 1)

    # Get some stuff
    self.s0 = self.beam.get_s0()
    self.m2 = self.gonio.get_rotation_axis()
    self.image_size = self.detector[0].get_image_size()

    # Get a random s1/phi
    i = uniform(0, self.image_size[0])
    j = uniform(1, self.image_size[1])
    self.s1 = matrix.col(self.detector[0].get_pixel_lab_coord((i, j)))
    self.s1 = self.s1.normalize() * matrix.col(self.s0).length()
    self.phi = uniform(0, 5)
    self.x0 = int(floor(i - 10))
    self.y0 = int(floor(j - 10))

    # Set some parameters
    self.sigma_divergence = self.beam.get_sigma_divergence(deg=False)
    self.delta_divergence = 3 * self.sigma_divergence
    self.grid_half_size = 4
    self.step_size = (self.delta_divergence / self.grid_half_size,
                      self.delta_divergence / self.grid_half_size)

    # Create the coordinate system
    self.cs = CoordinateSystem(self.m2, self.s0, self.s1, self.phi)

    # Create the map of s1 coordinates
    self.s1_map = transform.beam_vector_map(self.detector[0], self.beam, True)

    # Create the grid index generator
    self.generate_indices = transform.GridIndexGenerator(
        self.cs,
        self.x0,
        self.y0,
        self.step_size,
        self.grid_half_size,
        self.s1_map)

  def run(self):
    from scitbx import matrix

    for j in range(0, 20):
      for i in range(0, 20):

        xx = self.x0 + i
        yy = self.y0 + j
        if xx < 0 or yy < 0 or xx >= self.image_size[0] or yy >= self.image_size[0]:
          continue

        # Get the grid indices
        gi_1, gj_1 = self.generate_indices(j, i)

        # Get the grid indices
        xyz = matrix.col(self.detector[0].get_pixel_lab_coord(
            (self.x0 + i, self.y0 + j)))
        xyz = xyz.normalize() * matrix.col(self.s0).length()
        c1, c2 = matrix.col(self.cs.from_beam_vector(xyz))
        gi_2 = self.grid_half_size + c1 / self.step_size[0] + 0.5
        gj_2 = self.grid_half_size + c2 / self.step_size[1] + 0.5

        # Check both are the same
        eps = 1e-7
        try:
          assert(abs(gj_1 - gj_2) <= eps)
          assert(abs(gi_1 - gi_2) <= eps)
        except Exception, e:
          print gi_1, gi_2, gj_1, gj_2
          raise e

    # Test passed
    print 'OK'


if __name__ == '__main__':
  from dials.test import cd_auto
  with cd_auto(__file__):
    test = Test()
    test.run()
