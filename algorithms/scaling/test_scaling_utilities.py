'''
Tests for scaling utilities module.
'''
from math import sqrt, pi
import pytest
from mock import Mock
from libtbx.test_utils import approx_equal
from cctbx.sgtbx import space_group
from dials.array_family import flex
from dials.algorithms.scaling.scaling_utilities import \
  calc_crystal_frame_vectors, calc_theta_phi, create_sph_harm_table,\
  sph_harm_table, align_rotation_axis_along_z, parse_multiple_datasets,\
  set_wilson_outliers

@pytest.fixture(scope='module')
def test_space_group():
  """Create a space group object."""
  return space_group("C 2y")

@pytest.fixture(scope='module')
def mock_exp():
  """Create a mock experiments object."""
  exp = Mock()
  exp.beam.get_s0.return_value = (1.0, 0.0, 0.0)
  exp.goniometer.get_rotation_axis.return_value = (0.0, 0.0, 1.0)
  return exp

@pytest.fixture(scope='module')
def test_reflection_table():
  """Return a test reflection table."""
  return generate_reflection_table()

@pytest.fixture(scope='module')
def wilson_test_reflection_table():
  """Return a test reflection table."""
  rt = flex.reflection_table()
  rt['centric_flag'] = flex.bool([True, True, False, False])
  rt['Esq'] = flex.double([50.0, 10.0, 50.0, 10.0])
  return rt

def generate_reflection_table():
  """Create a reflection table with s1 and phi."""
  rt = flex.reflection_table()
  s1_vec = (1.0/sqrt(2.0), 0.0, 1.0/sqrt(2.0))
  rt['s1'] = flex.vec3_double([s1_vec, s1_vec, s1_vec])
  rt['phi'] = flex.double([0.0, 45.0, 90.0])
  return rt

def test_calc_crystal_frame_vectors(test_reflection_table, mock_exp):
  """Test the namesake function, to check that the vectors are correctly rotated
  into the crystal frame."""
  rt, exp = test_reflection_table, mock_exp
  s0_vec = (1.0, 0.0, 0.0)
  s1_vec = (1.0/sqrt(2.0), 0.0, 1.0/sqrt(2.0))
  reflection_table = calc_crystal_frame_vectors(rt, exp)
  assert list(reflection_table['s0']) == list(flex.vec3_double(
    [s0_vec, s0_vec, s0_vec]))
  assert approx_equal(list(reflection_table['s0c']), list(flex.vec3_double([
    s0_vec, (1.0/sqrt(2.0), -1.0/sqrt(2.0), 0.0), (0.0, -1.0, 0.0)])))
  assert approx_equal(list(reflection_table['s1c']), list(flex.vec3_double([
    s1_vec, (1.0/2.0, -1.0/2.0, 1.0/sqrt(2.0)),
    (0.0, -1.0/sqrt(2.0), 1.0/sqrt(2.0))])))

def test_align_rotation_axis_along_z():
  """Test the function to rotate the coordinate system such that the rotation
  axis is along z. In the test, the rotation axis is x, so we expect the
  transformation to be: x > z, y > y, z > -x, x+z > -x+z."""
  rot_axis = flex.vec3_double([(1.0, 0.0, 0.0)])
  vectors = flex.vec3_double([(1.0, 0.0, 0.0), (0.0, 1.0, 0.0), (0.0, 0.0, 1.0),
    (1.0, 0.0, 1.0)])
  rotated_vectors = align_rotation_axis_along_z(rot_axis, vectors)
  assert approx_equal(list(rotated_vectors), list(flex.vec3_double([
    (0.0, 0.0, 1.0), (0.0, 1.0, 0.0), (-1.0, 0.0, 0.0), (-1.0, 0.0, 1.0)])))

def test_sph_harm_table(test_reflection_table, mock_exp):
  """Simple test for the spherical harmonic table, constructing the table step
  by step, and verifying the values of a few easy-to-calculate entries.
  This also acts as a test for the calc_theta_phi function as well."""
  from scitbx import sparse # Needed to be able to assign to sph_h_t
  rt, exp = test_reflection_table, mock_exp
  reflection_table = calc_crystal_frame_vectors(rt, exp)
  theta_phi = calc_theta_phi(reflection_table['s0c'])
  assert approx_equal(list(theta_phi),
    [(pi/2.0, 0.0), (pi/2.0, -1.0*pi/4.0), (pi/2.0, -1.0*pi/2.0)])
  theta_phi_2 = calc_theta_phi(reflection_table['s1c'])
  assert approx_equal(list(theta_phi_2),
    [(pi/4.0, 0.0), (pi/4.0, -1.0*pi/4.0), (pi/4.0, -1.0*pi/2.0)])
  sph_h_t = create_sph_harm_table(theta_phi, theta_phi_2, 2)
  Y10 = ((3.0/(8.0*pi))**0.5)/2.0
  Y20 = -1.0*((5.0/(256.0*pi))**0.5)
  assert approx_equal(sph_h_t[0, 1], Y10)
  assert approx_equal(sph_h_t[1, 1], Y10)
  assert approx_equal(sph_h_t[2, 1], Y10)
  assert approx_equal(sph_h_t[0, 5], Y20)
  assert approx_equal(sph_h_t[1, 5], Y20)
  assert approx_equal(sph_h_t[2, 5], Y20)
  # Now test that you get the same by just calling the function.
  sht = sph_harm_table(rt, exp, 2)
  assert approx_equal(sht[0, 1], Y10)
  assert approx_equal(sht[1, 1], Y10)
  assert approx_equal(sht[2, 1], Y10)
  assert approx_equal(sht[0, 5], Y20)
  assert approx_equal(sht[1, 5], Y20)
  assert approx_equal(sht[2, 5], Y20)

def test_parse_multiple_datasets():
  """Test the namesake function. This expects a list of reflection tables, and
  selects on the column named 'id'."""
  rt1 = flex.reflection_table()
  rt1['id'] = flex.int([0, 0, 0])
  rt2 = flex.reflection_table()
  rt2['id'] = flex.int([1, 1, 2, 2])
  single_tables = parse_multiple_datasets([rt2])
  assert len(single_tables) == 2
  single_tables = parse_multiple_datasets([rt1, rt2])
  assert len(single_tables) == 3
  single_tables = parse_multiple_datasets([rt1])
  assert len(single_tables) == 1

def set_calculate_wilson_outliers(wilson_test_reflection_table):
  """Test the set wilson outliers function."""
  reflection_table = set_wilson_outliers(wilson_test_reflection_table)

  assert list(reflection_table.get_flags(
    reflection_table.flags.outlier_in_scaling)) == [True, False, True, False]
