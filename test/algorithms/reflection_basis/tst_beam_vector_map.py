from __future__ import absolute_import, division
from dials.algorithms.profile_model.gaussian_rs.transform import beam_vector_map

class Test(object):

  def __init__(self):
    from dials.model.serialize import load

    import libtbx.load_env
    try:
      dials_regression = libtbx.env.dist_path('dials_regression')
    except KeyError:
      print 'FAIL: dials_regression not configured'
      exit(0)

    import os

    filename = os.path.join(dials_regression,
        'centroid_test_data', 'sweep.json')

    # Load the sweep
    self.sweep = load.sweep(filename)

    # Get the models
    self.beam = self.sweep.get_beam()
    self.detector = self.sweep.get_detector()
    self.gonio = self.sweep.get_goniometer()
    self.scan = self.sweep.get_scan()

  def run(self):
    self.tst_at_corners()
    self.tst_sub_division_at_corners()
    self.tst_sub_division_at_centres()

  def tst_at_corners(self):

    from scitbx import matrix
    from random import randint
    assert(len(self.detector) == 1)

    # The detector beam vectors
    ds1 = beam_vector_map(self.detector[0], self.beam, True)
    expected_size = self.detector[0].get_image_size()[::-1]
    expected_size = tuple([e + 1 for e in expected_size])
    assert(ds1.all() == expected_size)

    s0 = self.beam.get_s0()
    m2 = self.gonio.get_rotation_axis()
    s0_length = matrix.col(self.beam.get_s0()).length()

    # Ensure a few random points are correct
    eps = 1e-7
    for k in range(1000):
      j = randint(0, ds1.all()[0]-1)
      i = randint(0, ds1.all()[1]-1)
      y = float(j)
      x = float(i)
      xyz = self.detector[0].get_pixel_lab_coord((x, y))
      s11 = matrix.col(xyz).normalize() * s0_length
      s12 = matrix.col(ds1[j,i])
      assert((s11 - s12).length() <= eps)

    # Test passed

  def tst_sub_division_at_corners(self):

    from scitbx import matrix
    from random import randint

    self.n_div = 2

    # The detector beam vectors
    ds1 = beam_vector_map(self.detector[0], self.beam, self.n_div, True)
    expected_size = self.detector[0].get_image_size()[::-1]
    expected_size = tuple([e * self.n_div + 1 for e in expected_size])
    assert(ds1.all() == expected_size)

    s0 = self.beam.get_s0()
    m2 = self.gonio.get_rotation_axis()
    s0_length = matrix.col(self.beam.get_s0()).length()

    # Ensure a few random points are correct
    eps = 1e-7
    for k in range(1000):
      j = randint(0, ds1.all()[0]-1)
      i = randint(0, ds1.all()[1]-1)
      y = float(j) / self.n_div
      x = float(i) / self.n_div
      xyz = self.detector[0].get_pixel_lab_coord((x, y))
      s11 = matrix.col(xyz).normalize() * s0_length
      s12 = matrix.col(ds1[j,i])
      assert((s11 - s12).length() <= eps)

    # Test passed

  def tst_sub_division_at_centres(self):

    from scitbx import matrix
    from random import randint

    self.n_div = 2

    # The detector beam vectors
    ds1 = beam_vector_map(self.detector[0], self.beam, self.n_div, False)
    expected_size = self.detector[0].get_image_size()[::-1]
    expected_size = tuple([e * self.n_div for e in expected_size])
    assert(ds1.all() == expected_size)

    s0 = self.beam.get_s0()
    m2 = self.gonio.get_rotation_axis()
    s0_length = matrix.col(self.beam.get_s0()).length()

    # Ensure a few random points are correct
    eps = 1e-7
    for k in range(1000):
      j = randint(0, ds1.all()[0]-1)
      i = randint(0, ds1.all()[1]-1)
      y = float(j + 0.5) / self.n_div
      x = float(i + 0.5) / self.n_div
      xyz = self.detector[0].get_pixel_lab_coord((x, y))
      s11 = matrix.col(xyz).normalize() * s0_length
      s12 = matrix.col(ds1[j,i])
      assert((s11 - s12).length() <= eps)

    # Test passed

if __name__ == '__main__':
  from dials.test import cd_auto
  with cd_auto(__file__):
    test = Test()
    test.run()
