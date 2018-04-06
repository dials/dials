from __future__ import absolute_import, division, print_function

import math
import random

from dials.algorithms.profile_model.gaussian_rs import *
import pytest

### Test the XDS coordinate system class

@pytest.fixture
def xdscoordinates():
  """Initialise the coordinate system"""
  # Coordinate sensitive to length of vectors, so need to ensure that
  # lengths of both s0 and s1 are equal
  coords = {
    's0': (0.013141995425357206,
           0.002199999234194632,
           1.4504754950989514),
    's1': (-0.01752795848400313,
           -0.24786554213968193,
           1.4290948735525306),
    'm2': (0.999975, -0.001289, -0.006968),
    'phi': 5.83575672475 * math.pi / 180
  }
  coords['cs'] = CoordinateSystem(coords['m2'], coords['s0'], coords['s1'], coords['phi'])
  return coords

def test_coordinate_system_data(xdscoordinates):
  """ Test all the input data """
  from scitbx import matrix

  eps = 1e-7
  s0 = matrix.col(xdscoordinates['s0'])
  s1 = matrix.col(xdscoordinates['s1'])
  m2 = matrix.col(xdscoordinates['m2'])
  assert(abs(matrix.col(xdscoordinates['cs'].s0()) - s0) <= eps)
  assert(abs(matrix.col(xdscoordinates['cs'].s0()) - s0) <= eps)
  assert(abs(matrix.col(xdscoordinates['cs'].s1()) - s1) <= eps)
  assert(abs(matrix.col(xdscoordinates['cs'].m2()) - m2.normalize()) <= eps)
  assert(abs(xdscoordinates['cs'].phi() - xdscoordinates['phi']) <= eps)

def test_axis_length(xdscoordinates):
  """Ensure axes are of unit length"""
  from scitbx import matrix

  # Check all axes are of length 1
  eps = 1e-7
  assert(abs(matrix.col(xdscoordinates['cs'].e1_axis()).length() - 1.0) <= eps)
  assert(abs(matrix.col(xdscoordinates['cs'].e2_axis()).length() - 1.0) <= eps)
  assert(abs(matrix.col(xdscoordinates['cs'].e3_axis()).length() - 1.0) <= eps)

def test_axis_orthogonal(xdscoordinates):
  """Ensure that the following are true:

  e1.s0 = 0, e1.s1 = 0
  e2.s1 = 0, e2.e1 = 0
  e3.e1 = 0, e3.p* = 0
  """
  from scitbx import matrix

  # Get as matrices
  e1 = matrix.col(xdscoordinates['cs'].e1_axis())
  e2 = matrix.col(xdscoordinates['cs'].e2_axis())
  e3 = matrix.col(xdscoordinates['cs'].e3_axis())
  s0 = matrix.col(xdscoordinates['s0'])
  s1 = matrix.col(xdscoordinates['s1'])

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

def test_limits(xdscoordinates):
  """ Test the coordinate system limits

  Ensure limits e1/e2 == |s1| and limit e3 == |s0 - s1|
  """
  from scitbx import matrix

  # Get the incident and diffracted beam vectors
  s0 = matrix.col(xdscoordinates['s0'])
  s1 = matrix.col(xdscoordinates['s1'])

  eps = 1e-7

  # Get the limits
  lim = xdscoordinates['cs'].limits()

  # Check the limits
  assert(abs(-1.0 - lim[0]) <= eps)
  assert(abs(1.0 - lim[1]) <= eps)


### Test the FromBeamVectorToXds class

@pytest.fixture
def beamvector():
  """Initialise the transform"""

  bv = {"s0": (0.013141995425357206,
                0.002199999234194632,
                1.4504754950989514),
        "s1": (-0.01752795848400313,
               -0.24786554213968193,
                1.4290948735525306),
        "m2": (0.999975, -0.001289, -0.006968),
        "phi": 5.83575672475 * math.pi / 180,
  }
  bv["cs"] = CoordinateSystem(bv["m2"], bv["s0"], bv["s1"], bv["phi"])
  return bv

def test_beamvector_coordinate_of_s1(beamvector):
  """Ensure that the coordinate of s1 is (0, 0, 0)"""

  # Get the coordinate at s1
  s_dash = beamvector["s1"]
  c1, c2 = beamvector["cs"].from_beam_vector(s_dash)

  # Ensure that it is at the origin
  assert c1 == pytest.approx(0.0)
  assert c2 == pytest.approx(0.0)

def test_beamvector_limit(beamvector):
  """ Calculate the coordinate at the limits.

  Ensure that coordinate where s1' is orthogonal to s1 is at limit.
  """
  from scitbx import matrix

  # Get the limit of s1'
  s_dash = matrix.col(beamvector["s1"]).cross(matrix.col(beamvector["s0"]))
  s_dash = s_dash.normalize() * matrix.col(beamvector["s1"]).length()

  # Rotate arbitrarily
  s_dash = s_dash.rotate_around_origin(matrix.col(beamvector["s1"]), random.uniform(0, 360), deg=True)

  # Get the c1, c2 coordinate
  c1, c2 = beamvector["cs"].from_beam_vector(s_dash)

  # Check the point is equal to the limit in rs
  assert math.sqrt(c1**2 + c2**2) == pytest.approx(abs(beamvector["cs"].limits()[0]))


### Test the TestFromRotationAngle class

@pytest.fixture
def rotationangle():
  """Initialise the transform"""

  ra = {
      "s0": (0.013141995425357206,
              0.002199999234194632,
              1.4504754950989514),
      "s1": (-0.01752795848400313,
             -0.24786554213968193,
              1.4290948735525306),
      "m2": (0.999975, -0.001289, -0.006968),
      "phi": 5.83575672475 * math.pi / 180,
  }

  ra["cs"] = CoordinateSystem(ra["m2"], ra["s0"], ra["s1"], ra["phi"])
  return ra

def test_from_rotation_angle_coordinate_of_phi(rotationangle):
  """Ensure that the coordinate of s1 is (0, 0, 0)"""

  # Get the coordinate at phi
  phi_dash = rotationangle['phi']
  c3 = rotationangle["cs"].from_rotation_angle(phi_dash)

  # Ensure that it is at the origin
  assert c3 == pytest.approx(0.0)

def test_from_rotation_angle_e3_coordinate_approximation(rotationangle):
  # Select a random rotation from phi
  s_dash = rotationangle["s1"]
  phi_dash = rotationangle["phi"] + (2.0*random.random() - 1.0) * math.pi / 180

  # Calculate the XDS coordinate, this class uses an approximation
  # for c3 = zeta * (phi' - phi)
  c3 = rotationangle["cs"].from_rotation_angle(phi_dash)
  c3_2 = rotationangle["cs"].from_rotation_angle_fast(phi_dash)

  # Check the true and approximate value are almost equal to 4dp
  assert c3 == pytest.approx(c3_2, abs=1e-4)


### Test the ToBeamVector class

def test_to_beamvector_xds_origin(beamvector):
  """Test the beam vector at the XDS origin is equal to s1."""
  from scitbx import matrix
  eps = 1e-7
  s_dash = beamvector["cs"].to_beam_vector((0, 0))
  assert(abs(matrix.col(s_dash) - matrix.col(beamvector["s1"])) <= eps)

def test_to_beamvector_far_out_coordinates(beamvector):
  """Test some large coordinates, 1 valid and the other invalid. (i.e.
  a coordinate that cannot be mapped onto the ewald sphere)."""
  eps = 1e-7

  # Setting c2 and c3 to zero
  c2 = 0

  # A large value which is still valid
  c1 = 1.0 - eps
  s_dash = beamvector["cs"].to_beam_vector((c1, c2))

  # A large value which is raises an exception
  with pytest.raises(RuntimeError):
    c1 = 1.0 + eps
    s_dash = beamvector["cs"].to_beam_vector((c1, c2))

  # Setting c2 and c3 to zero
  c1 = 0
  c3 = 0

  # A large value which is still valid
  c2 = 1.0 - eps
  s_dash = beamvector["cs"].to_beam_vector((c1, c2))

  # A large value which is raises an exception
  with pytest.raises(RuntimeError):
    c2 = 1.0 + eps
    s_dash = beamvector["cs"].to_beam_vector((c1, c2))

def test_to_beamvector_forward_and_reverse_transform(beamvector):
  """Test the forward and reverse Beam Vector -> XDS transforms Create
  a beam vector, transform it to XDS and then transform back. The new
  value should be equal to the original value."""

  from scitbx import matrix

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
    s_dash = matrix.col(beamvector["s1"]) + matrix.col((random_shift(),
                                               random_shift(),
                                               random_shift()))
    s_dash = s_dash.normalize() * matrix.col(beamvector["s1"]).length()

    # Calculate the XDS coordinate of the vector
    c1, c2 = beamvector["cs"].from_beam_vector(s_dash)

    # Calculate the beam vector from the XDS coordinate
    s_dash_2 = beamvector["cs"].to_beam_vector((c1, c2))

    # Check the vectors are almost equal
    assert(abs(matrix.col(s_dash) - matrix.col(s_dash_2)) <= eps)

### Test the TestToRotationAngle class

def test_forward_and_backward(rotationangle):
  # Set the parameters
  min_shift = -20.0 * math.pi / 180.0
  max_shift = +20.0 * math.pi / 180.0
  range_shift = max_shift - min_shift
  def random_shift():
    return min_shift + random.random() * range_shift

  # Loop a number of times
  num = 1000
  for i in range(num):

    # Create a rotation angle
    phi_dash = rotationangle["phi"] + random_shift()

    # Calculate the XDS coordinate of the vector
    c3 = rotationangle["cs"].from_rotation_angle(phi_dash)

    # Calculate the beam vector from the XDS coordinate
    phi_dash_2 = rotationangle["cs"].to_rotation_angle(c3)

    # Check the vectors are almost equal
    assert phi_dash_2 == pytest.approx(phi_dash)

def test_origin(rotationangle):
  phi_dash = rotationangle["cs"].to_rotation_angle(0.0)
  assert phi_dash == pytest.approx(rotationangle["phi"])

def test_far_out_coordinates(rotationangle):
  """Test some large coordinates, 1 valid and the other invalid. (i.e.
  a coordinate that cannot be mapped onto the ewald sphere)."""

  eps = 1e-7

  # Get the limits
  lim = rotationangle["cs"].limits()[2:]
  lim0 = max(lim)
  lim1 = min(lim)

  # Setting c2 and c3 to zero
  c3 = lim0 - eps

  # A large value which is still valid
  phi_dash = rotationangle["cs"].to_rotation_angle(c3)

  # Setting c2 and c3 to zero
  c3 = lim1 + eps

  # A large value which is still valid
  phi_dash = rotationangle["cs"].to_rotation_angle(c3)

  # A large value which is raises an exception
  with pytest.raises(RuntimeError):
    c3 = lim0 + eps
    phi_dash = rotationangle["cs"].to_rotation_angle(c3)
    print(phi_dash)

  with pytest.raises(RuntimeError):
    c3 = lim1 - eps
    phi_dash = rotationangle["cs"].to_rotation_angle(c3)

def test_e3_coordinate_approximation(rotationangle):
  # Select a random rotation from phi
  phi_dash = rotationangle["phi"] + (2.0*random.random() - 1.0) * math.pi / 180

  # Calculate the XDS coordinate, this class uses an approximation
  # for c3 = zeta * (phi' - phi)
  c3 = rotationangle["cs"].from_rotation_angle(phi_dash)
  phi_dash_2 = rotationangle["cs"].to_rotation_angle_fast(c3)

  # Check the true and approximate value are almost equal to 4dp
  assert phi_dash == pytest.approx(phi_dash_2, abs=1e-4)
