from __future__ import absolute_import, division
from dials.algorithms.profile_model.gaussian_rs import *

class TestCoordinateSystem(object):
  """Test the XDS coordinate system class"""

  def __init__(self):
    """Initialise the coordinate system"""

    from math import pi

    # Coordinate sensitive to length of vectors, so need to ensure that
    # lengths of both s0 and s1 are equal
    self.s0 = (0.013141995425357206,
                0.002199999234194632,
                1.4504754950989514)
    self.s1 = (-0.01752795848400313,
               -0.24786554213968193,
                1.4290948735525306)
    self.m2 = (0.999975, -0.001289, -0.006968)
    self.phi = 5.83575672475 * pi / 180
    self.cs = CoordinateSystem(self.m2, self.s0, self.s1, self.phi)

  def test_data(self):
    """ Test all the input data """

    from scitbx import matrix

    eps = 1e-7
    s0 = matrix.col(self.s0)
    s1 = matrix.col(self.s1)
    m2 = matrix.col(self.m2)
    assert(abs(matrix.col(self.cs.s0()) - s0) <= eps)
    assert(abs(matrix.col(self.cs.s1()) - s1) <= eps)
    assert(abs(matrix.col(self.cs.m2()) - m2.normalize()) <= eps)
    assert(abs(self.cs.phi() - self.phi) <= eps)
    print 'OK'

  def test_axis_length(self):
    """Ensure axes are of unit length"""
    from scitbx import matrix

    # Check all axes are of length 1
    eps = 1e-7
    assert(abs(matrix.col(self.cs.e1_axis()).length() - 1.0) <= eps)
    assert(abs(matrix.col(self.cs.e2_axis()).length() - 1.0) <= eps)
    assert(abs(matrix.col(self.cs.e3_axis()).length() - 1.0) <= eps)

    # Test passed
    print "OK"

  def test_axis_orthogonal(self):
    """Ensure that the following are true:

    e1.s0 = 0, e1.s1 = 0
    e2.s1 = 0, e2.e1 = 0
    e3.e1 = 0, e3.p* = 0

    """
    from scitbx import matrix

    # Get as matrices
    e1 = matrix.col(self.cs.e1_axis())
    e2 = matrix.col(self.cs.e2_axis())
    e3 = matrix.col(self.cs.e3_axis())
    s0 = matrix.col(self.s0)
    s1 = matrix.col(self.s1)

    eps = 1e-7

    # Check that e1 is orthogonal to s0 and s1
    assert(abs(e1.dot(s0) - 0.0) <= eps)
    assert(abs(e1.dot(s1) - 0.0) <= eps)

    # Check that e2 is orthogonal to s1 and e1
    assert(abs(e2.dot(s1) - 0.0) <= eps)
    assert(abs(e2.dot(e1) - 0.0) <= eps)

    # Check that e3 is orthogonal to e1 and p* = s1 - s0
    assert(abs(e3.dot(e1) - 0.0) <= eps)
    assert(abs(e3.dot(s1 - s0) - 0.0) <= eps)

    # Test passed
    print "OK"

  def test_limits(self):
    """ Test the coordinate system limits

    Ensure limits e1/e2 == |s1| and limit e3 == |s0 - s1|

    """
    from scitbx import matrix

    # Get the incident and diffracted beam vectors
    s0 = matrix.col(self.s0)
    s1 = matrix.col(self.s1)

    eps = 1e-7

    # Get the limits
    lim = self.cs.limits()

    # Check the limits
    assert(abs(-1.0 - lim[0]) <= eps)
    assert(abs(1.0 - lim[1]) <= eps)

    # Test passed
    print 'OK'

  def __call__(self):
    """Run the tests"""
    self.test_data()
    self.test_axis_length()
    self.test_axis_orthogonal()
    self.test_limits()


class TestFromBeamVector(object):
  """Test the FromBeamVectorToXds class"""

  def __init__(self):
    """Initialise the transform"""

    from math import pi

    self.s0 = (0.013141995425357206,
                0.002199999234194632,
                1.4504754950989514)
    self.s1 = (-0.01752795848400313,
               -0.24786554213968193,
                1.4290948735525306)
    self.m2 = (0.999975, -0.001289, -0.006968)
    self.phi = 5.83575672475 * pi / 180

    self.cs = CoordinateSystem(self.m2, self.s0, self.s1, self.phi)

  def test_coordinate_of_s1(self):
    """Ensure that the coordinate of s1 is (0, 0, 0)"""

    eps = 1e-7

    # Get the coordinate at s1
    s_dash = self.s1
    c1, c2 = self.cs.from_beam_vector(s_dash)

    # Ensure that it is at the origin
    assert(abs(c1 - 0.0) <= eps)
    assert(abs(c2 - 0.0) <= eps)

    # Test passed
    print "OK"

  def test_limit(self):
    """ Calculate the coordinate at the limits.

    Ensure that coordinate where s1' is orthogonal to s1 is at limit.

    """
    from scitbx import matrix
    from math import sqrt
    from random import uniform
    eps = 1e-7

    # Get the limit of s1'
    s_dash = matrix.col(self.s1).cross(matrix.col(self.s0))
    s_dash = s_dash.normalize() * matrix.col(self.s1).length()

    # Rotate arbitrarily
    s_dash = s_dash.rotate(matrix.col(self.s1), uniform(0, 360), deg=True)

    # Get the c1, c2 coordinate
    c1, c2 = self.cs.from_beam_vector(s_dash)

    # Check the point is equal to the limit in rs
    assert(abs(sqrt(c1**2 + c2**2) - abs(self.cs.limits()[0])) <= eps)

    # Test passed
    print 'OK'

  def __call__(self):
    """Run all the tests"""
    self.test_coordinate_of_s1()
    self.test_limit()


class TestFromRotationAngle(object):
  """Test the TestFromRotationAngle class"""

  def __init__(self):
    """Initialise the transform"""

    from math import pi

    self.s0 = (0.013141995425357206,
                0.002199999234194632,
                1.4504754950989514)
    self.s1 = (-0.01752795848400313,
               -0.24786554213968193,
                1.4290948735525306)
    self.m2 = (0.999975, -0.001289, -0.006968)
    self.phi = 5.83575672475 * pi / 180


    self.cs = CoordinateSystem(self.m2, self.s0, self.s1, self.phi)

    self.p_star = self.cs.p_star()

  def test_coordinate_of_phi(self):
    """Ensure that the coordinate of s1 is (0, 0, 0)"""

    eps = 1e-7

    # Get the coordinate at phi
    phi_dash = self.phi
    c3 = self.cs.from_rotation_angle(phi_dash)

    # Ensure that it is at the origin
    assert(abs(c3 - 0.0) <= eps)

    # Test passed
    print "OK"

  def test_e3_coordinate_approximation(self):

    from math import pi
    import random

    eps = 1e-4

    # Select a random rotation from phi
    s_dash = self.s1
    phi_dash = self.phi + (2.0*random.random() - 1.0) * pi / 180

    # Calculate the XDS coordinate, this class uses an approximation
    # for c3 = zeta * (phi' - phi)
    c3 = self.cs.from_rotation_angle(phi_dash)
    c3_2 = self.cs.from_rotation_angle_fast(phi_dash)

    # Check the true and approximate value are almost equal to 4dp
    assert(abs(c3 - c3_2) < eps)

    # Test passed
    print "OK"

  def __call__(self):
    """Run all the tests"""
    self.test_coordinate_of_phi()
    self.test_e3_coordinate_approximation()


class TestToBeamVector(object):
  """Test the ToBeamVector class"""

  def __init__(self):
    """Initialise the transform"""

    from math import pi

    self.s0 = (0.013141995425357206,
                0.002199999234194632,
                1.4504754950989514)
    self.s1 = (-0.01752795848400313,
               -0.24786554213968193,
                1.4290948735525306)
    self.m2 = (0.999975, -0.001289, -0.006968)
    self.phi = 5.83575672475 * pi / 180

    self.cs = CoordinateSystem(self.m2, self.s0, self.s1, self.phi)

  def test_xds_origin(self):
    """Test the beam vector at the XDS origin is equal to s1."""
    from scitbx import matrix
    eps = 1e-7
    s_dash = self.cs.to_beam_vector((0, 0))
    assert(abs(matrix.col(s_dash) - matrix.col(self.s1)) <= eps)

    print "OK"

  def test_far_out_coordinates(self):
    """Test some large coordinates, 1 valid and the other invalid. (i.e.
    a coordinate that cannot be mapped onto the ewald sphere)."""


    eps = 1e-7

    # Setting c2 and c3 to zero
    c2 = 0

    # A large value which is still valid
    c1 = 1.0 - eps
    s_dash = self.cs.to_beam_vector((c1, c2))

    # A large value which is raises an exception
    try:
      c1 = 1.0 + eps
      s_dash = self.cs.to_beam_vector((c1, c2))
      assert(False)
    except RuntimeError:

      #Test passed
      print "OK"

    # Setting c2 and c3 to zero
    c1 = 0
    c3 = 0

    # A large value which is still valid
    c2 = 1.0 - eps
    s_dash = self.cs.to_beam_vector((c1, c2))

    # A large value which is raises an exception
    try:
      c2 = 1.0 + eps
      s_dash = self.cs.to_beam_vector((c1, c2))
      assert(False)
    except RuntimeError:

      #Test passed
      print "OK"

  def test_forward_and_reverse_transform(self):
    """Test the forward and reverse Beam Vector -> XDS transforms Create
    a beam vector, transform it to XDS and then transform back. The new
    value should be equal to the original value."""

    from scitbx import matrix
    import random

    eps = 1e-7

    # Set the parameters
    min_shift = -0.5
    max_shift = +0.5
    range_shift = max_shift - min_shift
    def random_shift():
      return min_shift + random.random() * range_shift

    # Loop a number of times
    num = 1000
    for i in range(num):

      # Create a beam vector
      s_dash = matrix.col(self.s1) + matrix.col((random_shift(),
                                                 random_shift(),
                                                 random_shift()))
      s_dash = s_dash.normalize() * matrix.col(self.s1).length()

      # Calculate the XDS coordinate of the vector
      c1, c2 = self.cs.from_beam_vector(s_dash)

      # Calculate the beam vector from the XDS coordinate
      s_dash_2 = self.cs.to_beam_vector((c1, c2))

      # Check the vectors are almost equal
      assert(abs(matrix.col(s_dash) - matrix.col(s_dash_2)) <= eps)

    # Test passed
    print "OK"

  def __call__(self):
    """Run the tests"""
    self.test_xds_origin()
    self.test_far_out_coordinates()
    self.test_forward_and_reverse_transform()


class TestToRotationAngle(object):
  """Test the TestToRotationAngle class"""

  def __init__(self):
    """Initialise the transform"""

    from math import pi

    self.s0 = (0.013141995425357206,
                0.002199999234194632,
                1.4504754950989514)
    self.s1 = (-0.01752795848400313,
               -0.24786554213968193,
                1.4290948735525306)
    self.m2 = (0.999975, -0.001289, -0.006968)
    self.phi = 5.83575672475 * pi / 180

    self.cs = CoordinateSystem(self.m2, self.s0, self.s1, self.phi)

  def test_forward_and_backward(self):

    from math import pi
    import random
    eps = 1e-7

    # Set the parameters
    min_shift = -20.0 * pi / 180.0
    max_shift = +20.0 * pi / 180.0
    range_shift = max_shift - min_shift
    def random_shift():
      return min_shift + random.random() * range_shift

    # Loop a number of times
    num = 1000
    for i in range(num):

      # Create a rotation angle
      phi_dash = self.phi + random_shift()

      # Calculate the XDS coordinate of the vector
      c3 = self.cs.from_rotation_angle(phi_dash)

      # Calculate the beam vector from the XDS coordinate
      phi_dash_2 = self.cs.to_rotation_angle(c3)

      # Check the vectors are almost equal
      assert(abs(phi_dash_2 - phi_dash) <= eps)

    # Test passed
    print "OK"

  def test_origin(self):
    eps = 1e-7
    phi_dash = self.cs.to_rotation_angle(0.0)
    assert(abs(phi_dash - self.phi) <= eps)
    print 'OK'

  def test_far_out_coordinates(self):
    """Test some large coordinates, 1 valid and the other invalid. (i.e.
    a coordinate that cannot be mapped onto the ewald sphere)."""

    eps = 1e-7

    # Get the limits
    lim = self.cs.limits()[2:]
    lim0 = max(lim)
    lim1 = min(lim)

    # Setting c2 and c3 to zero
    c3 = lim0 - eps

    # A large value which is still valid
    phi_dash = self.cs.to_rotation_angle(c3)

    # Setting c2 and c3 to zero
    c3 = lim1 + eps

    # A large value which is still valid
    phi_dash = self.cs.to_rotation_angle(c3)

    # A large value which is raises an exception
    try:
      c3 = lim0 + eps
      phi_dash = self.cs.to_rotation_angle(c3)
      print phi_dash
      assert(False)
    except RuntimeError:
      pass

    try:
      c3 = lim1 - eps
      phi_dash = self.cs.to_rotation_angle(c3)
      assert(False)
    except RuntimeError:
      pass

    #Test passed
    print "OK"

  def test_e3_coordinate_approximation(self):

    from math import pi
    import random

    eps = 1e-4

    # Select a random rotation from phi
    phi_dash = self.phi + (2.0*random.random() - 1.0) * pi / 180

    # Calculate the XDS coordinate, this class uses an approximation
    # for c3 = zeta * (phi' - phi)
    c3 = self.cs.from_rotation_angle(phi_dash)
    phi_dash_2 = self.cs.to_rotation_angle_fast(c3)

    # Check the true and approximate value are almost equal to 4dp
    assert(abs(phi_dash - phi_dash_2) < eps)

    # Test passed
    print "OK"

  def __call__(self):
    """Run all the tests"""
    self.test_forward_and_backward()
    self.test_origin()
    self.test_far_out_coordinates()
    self.test_e3_coordinate_approximation()


class Test(object):
  def __init__(self):
    self.tst_coordinate_system = TestCoordinateSystem()
    self.tst_from_beam_vector = TestFromBeamVector()
    self.tst_from_rotation_angle = TestFromRotationAngle()
    self.tst_to_beam_vector = TestToBeamVector()
    self.tst_to_rotation_angle = TestToRotationAngle()

  def run(self):
    self.tst_coordinate_system()
    self.tst_from_beam_vector()
    self.tst_from_rotation_angle()
    self.tst_to_beam_vector()
    self.tst_to_rotation_angle()


if __name__ == '__main__':
  from dials.test import cd_auto
  with cd_auto(__file__):
    test = Test()
    test.run()
