from __future__ import absolute_import, division

class TestMapFramesForward(object):

  def __init__(self, filename):

    from math import pi
    from dials.model.serialize import load
    from dials.algorithms.profile_model.gaussian_rs.transform import MapFramesForward
    from dials.algorithms.profile_model.gaussian_rs import BBoxCalculator3D

    # Load the sweep
    self.sweep = load.sweep(filename)

    # Get the models
    self.beam = self.sweep.get_beam()
    self.detector = self.sweep.get_detector()
    self.gonio = self.sweep.get_goniometer()
    self.scan = self.sweep.get_scan()

    # Set the delta_divergence/mosaicity
    self.n_sigma = 3
    self.sigma_divergence = 0.060 * pi / 180
    self.mosaicity = 0.154 * pi / 180
    self.delta_divergence = self.n_sigma * self.sigma_divergence
    self.delta_mosaicity = self.n_sigma * self.mosaicity

    # Set the grid size
    self.grid_size = (4, 4, 4)

    # Create the E3 fraction object
    self.transform = MapFramesForward(
        self.scan.get_array_range()[0],
        self.scan.get_oscillation(deg=False)[0],
        self.scan.get_oscillation(deg=False)[1],
        self.mosaicity,
        self.n_sigma,
        self.grid_size[2])

    # Create the bounding box calculator
    self.calculate_bbox = BBoxCalculator3D(
        self.beam, self.detector, self.gonio, self.scan,
        self.delta_divergence,
        self.delta_mosaicity)

  def __call__(self):

    from scitbx import matrix
    from random import uniform
    from dials.algorithms.profile_model.gaussian_rs import CoordinateSystem
    from scitbx.array_family import flex
    assert(len(self.detector) == 1)
    s0 = self.beam.get_s0()
    m2 = self.gonio.get_rotation_axis()
    s0_length = matrix.col(self.beam.get_s0()).length()

    for i in range(100):

      # Get random x, y, z
      x = uniform(0, 2000)
      y = uniform(0, 2000)
      z = uniform(0, 9)

      # Get random s1, phi, panel
      s1 = matrix.col(self.detector[0].get_pixel_lab_coord(
          (x, y))).normalize() * s0_length
      phi = self.scan.get_angle_from_array_index(z, deg=False)
      panel = 0

      # Calculate the bounding box
      bbox = self.calculate_bbox(s1, z, panel)
      x1, x2 = bbox[0], bbox[1]
      y1, y2 = bbox[2], bbox[3]
      z1, z2 = bbox[4], bbox[5]

      # Create the XDS coordinate system
      xcs = CoordinateSystem(m2, s0, s1, phi)

      # Calculate the transform fraction
      fraction = self.transform(bbox[4:], phi, xcs.zeta())

      # Ensure the minimum and maximum are 0 < 1
      fmax = flex.max(fraction)
      fmin = flex.min(fraction)
      assert(fmax <= (1.0 + 5e-15) and fmax > 0.0), "%.16f not between 0 and 1" % fmax
      assert(fmin >= 0.0 and fmin <= 1.0)

      # Ensure the fraction for each image frame adds up to 1.0 for
      # all those frames completely within the grid
      for j in range(1, fraction.all()[0]-1):
        tot = flex.sum(fraction[j:j+1,:])
        assert(abs(tot - 1.0) < 1e-7)

      # Ensure the frames follow a progression through the grid. I.e,
      # check that values increase then decrease and don't jump around
      for j in range(fraction.all()[0]):
        f = fraction[j:j+1,:]
        last = f[0]
        rev = False
        for i in range(1,len(f)):
          curr = f[1]
          if rev == False:
            if curr < last:
              rev = True
          else:
            assert(curr <= last)
          last = curr

    # Test passed


class TestMapFramesReverse(object):

  def __init__(self, filename):

    from math import pi
    from dials.model.serialize import load
    from dials.algorithms.profile_model.gaussian_rs.transform import MapFramesReverse
    from dials.algorithms.profile_model.gaussian_rs import BBoxCalculator3D

    # Load the sweep
    self.sweep = load.sweep(filename)

    # Get the models
    self.beam = self.sweep.get_beam()
    self.detector = self.sweep.get_detector()
    self.gonio = self.sweep.get_goniometer()
    self.scan = self.sweep.get_scan()

    # Set the delta_divergence/mosaicity
    self.n_sigma = 3
    self.sigma_divergence = 0.060 * pi / 180
    self.mosaicity = 0.154 * pi / 180
    self.delta_divergence = self.n_sigma * self.sigma_divergence
    self.delta_mosaicity = self.n_sigma * self.mosaicity

    # Set the grid size
    self.grid_size = (4, 4, 4)

    # Create the E3 fraction object
    self.transform = MapFramesReverse(
        self.scan.get_array_range()[0],
        self.scan.get_oscillation(deg=False)[0],
        self.scan.get_oscillation(deg=False)[1],
        self.mosaicity,
        self.n_sigma,
        self.grid_size[2])

    # Create the bounding box calculator
    self.calculate_bbox = BBoxCalculator3D(
        self.beam, self.detector, self.gonio, self.scan,
        self.delta_divergence,
        self.delta_mosaicity)

  def __call__(self):

    from scitbx import matrix
    from random import uniform
    from dials.algorithms.profile_model.gaussian_rs import CoordinateSystem
    from scitbx.array_family import flex

    s0 = self.beam.get_s0()
    m2 = self.gonio.get_rotation_axis()
    s0_length = matrix.col(self.beam.get_s0()).length()

    for i in range(100):

      # Get random x, y, z
      x = uniform(0, 2000)
      y = uniform(0, 2000)
      z = uniform(0, 9)

      # Get random s1, phi, panel
      s1 = matrix.col(self.detector[0].get_pixel_lab_coord(
          (x, y))).normalize() * s0_length
      phi = self.scan.get_angle_from_array_index(z, deg=False)
      panel = 0

      # Calculate the bounding box
      bbox = self.calculate_bbox(s1, phi, panel)
      x1, x2 = bbox[0], bbox[1]
      y1, y2 = bbox[2], bbox[3]
      z1, z2 = bbox[4], bbox[5]
      if x1 == 0 or y1 == 0 or z1 == 0:
        continue
      if x2 == 2000 or y2 == 2000 or z2 == 9:
        continue

      # Create the XDS coordinate system
      xcs = CoordinateSystem(m2, s0, s1, phi)

      # Calculate the transform fraction
      fraction = self.transform(bbox[4:], phi, xcs.zeta())

      # Ensure the minimum and maximum are 0 < 1
      fmax = flex.max(fraction)
      fmin = flex.min(fraction)
      assert(fmax <= 1.0 and fmax > 0.0)
      assert(fmin >= 0.0 and fmin <= 1.0)

      # Ensure the fraction for image adds up to 1.0 for
      # all those images completely within the image
      for v3 in range(fraction.all()[0]):
        tot = flex.sum(fraction[v3:v3+1,:])
        assert(abs(tot - 1.0) < 1e-7)

      # Ensure the frames follow a progression through the grid. I.e,
      # check that values increase then decrease and don't jump around
      for v3 in range(fraction.all()[0]):
        f = fraction[v3:v3+1,:]
        last = f[0]
        rev = False
        for i in range(1,len(f)):
          curr = f[1]
          if rev == False:
            if curr < last:
              rev = True
          else:
            assert(curr <= last)
          last = curr

    # Test passed


class TestMapForwardReverse(object):

  def __init__(self, filename):

    from math import pi
    from dials.model.serialize import load
    from dials.algorithms.profile_model.gaussian_rs.transform import MapFramesReverse
    from dials.algorithms.profile_model.gaussian_rs.transform import MapFramesForward
    from dials.algorithms.profile_model.gaussian_rs import BBoxCalculator3D

    # Load the sweep
    self.sweep = load.sweep(filename)

    # Get the models
    self.beam = self.sweep.get_beam()
    self.detector = self.sweep.get_detector()
    self.gonio = self.sweep.get_goniometer()
    self.scan = self.sweep.get_scan()

    # Set the delta_divergence/mosaicity
    self.n_sigma = 3
    self.sigma_divergence = 0.060 * pi / 180
    self.mosaicity = 0.154 * pi / 180
    self.delta_divergence = self.n_sigma * self.sigma_divergence
    self.delta_mosaicity = self.n_sigma * self.mosaicity

    # Set the grid size
    self.grid_size = (4, 4, 4)

    # Create the E3 fraction object
    self.transform_forward = MapFramesForward(
        self.scan.get_array_range()[0],
        self.scan.get_oscillation(deg=False)[0],
        self.scan.get_oscillation(deg=False)[1],
        self.mosaicity,
        self.n_sigma,
        self.grid_size[2])

    # Create the E3 fraction object
    self.transform_reverse = MapFramesReverse(
        self.scan.get_array_range()[0],
        self.scan.get_oscillation(deg=False)[0],
        self.scan.get_oscillation(deg=False)[1],
        self.mosaicity,
        self.n_sigma,
        self.grid_size[2])

    # Create the bounding box calculator
    self.calculate_bbox = BBoxCalculator3D(
        self.beam, self.detector, self.gonio, self.scan,
        self.delta_divergence,
        self.delta_mosaicity)

  def __call__(self):

    from scitbx import matrix
    from random import uniform
    from dials.algorithms.profile_model.gaussian_rs import CoordinateSystem

    s0 = self.beam.get_s0()
    m2 = self.gonio.get_rotation_axis()
    s0_length = matrix.col(self.beam.get_s0()).length()

    for i in range(100):

      # Get random x, y, z
      x = uniform(0, 2000)
      y = uniform(0, 2000)
      z = uniform(0, 9)

      # Get random s1, phi, panel
      s1 = matrix.col(self.detector[0].get_pixel_lab_coord(
          (x, y))).normalize() * s0_length
      phi = self.scan.get_angle_from_array_index(z, deg=False)
      panel = 0

      # Calculate the bounding box
      bbox = self.calculate_bbox(s1, phi, panel)
      x1, x2 = bbox[0], bbox[1]
      y1, y2 = bbox[2], bbox[3]
      z1, z2 = bbox[4], bbox[5]

      # Create the XDS coordinate system
      xcs = CoordinateSystem(m2, s0, s1, phi)

      # Calculate the transform fraction
      forward_fraction = self.transform_forward(bbox[4:], phi, xcs.zeta())

      # Calculate the transform fraction
      reverse_fraction = self.transform_reverse(bbox[4:], phi, xcs.zeta())

      # Check the same points are non-zero
      eps = 1e-7
      for j in range(forward_fraction.all()[0]):
        for i in range(forward_fraction.all()[1]):
          if forward_fraction[j,i] > 0.0:
            assert(reverse_fraction[i,j] > 0.0)
          else:
            assert(reverse_fraction[i,j] < eps)

    # Test passed


class Test(object):

  def __init__(self):

    import libtbx.load_env
    try:
      dials_regression = libtbx.env.dist_path('dials_regression')
    except KeyError:
      print 'FAIL: dials_regression not configured'
      exit(0)

    import os

    filename = os.path.join(dials_regression,
        'centroid_test_data', 'sweep.json')

    self.tst_map_frames_forward = TestMapFramesForward(filename)
    self.tst_map_frames_reverse = TestMapFramesReverse(filename)
    self.tst_map_forward_reverse = TestMapForwardReverse(filename)

  def run(self):

    self.tst_map_frames_forward()
    self.tst_map_frames_reverse()
    self.tst_map_forward_reverse()

if __name__ == '__main__':
  from dials.test import cd_auto
  with cd_auto(__file__):
    test = Test()
    test.run()
