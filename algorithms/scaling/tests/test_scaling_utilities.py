'''
Tests for scaling utilities module.
'''
from math import sqrt, pi
from dials.array_family import flex
from libtbx.test_utils import approx_equal
from dxtbx.model.experiment_list import ExperimentList
from dxtbx.model import Crystal, Scan, Beam, Goniometer, Detector, Experiment
from dials.algorithms.scaling.scaling_utilities import calc_crystal_frame_vectors,\
  calc_theta_phi, create_sph_harm_table, sph_harm_table, align_rotation_axis_along_z, \
  parse_multiple_datasets
#import pytest

def generated_single_exp():
  """Generate an experiment."""
  experiments = ExperimentList()
  exp_dict = {"__id__" : "crystal", "real_space_a": [1.0, 0.0, 0.0],
              "real_space_b": [0.0, 1.0, 0.0], "real_space_c": [0.0, 0.0, 2.0],
              "space_group_hall_symbol": " C 2y"}
  crystal = Crystal.from_dict(exp_dict)
  scan = Scan(image_range=[0, 90], oscillation=[0.0, 1.0])
  beam = Beam(s0=(1.0, 0.0, 0.0))
  goniometer = Goniometer((0.0, 0.0, 1.0))
  detector = Detector()
  experiments.append(Experiment(beam=beam, scan=scan, goniometer=goniometer,
    detector=detector, crystal=crystal))
  return experiments

def generate_input_vectors():
  """Generate the input for the following tests."""
  rt = flex.reflection_table()
  s1_vec = (1.0/sqrt(2.0), 0.0, 1.0/sqrt(2.0))
  rt['s1'] = flex.vec3_double([s1_vec, s1_vec, s1_vec])
  rt['phi'] = flex.double([0.0, 45.0, 90.0])
  exp = generated_single_exp()[0]
  return rt, exp

def test_calc_crystal_frame_vectors():
  """Test the namesake function, to check that the vectors are correctly rotated
  into the crystal frame."""
  rt, exp = generate_input_vectors()
  s0_vec = (1.0, 0.0, 0.0)
  s1_vec = (1.0/sqrt(2.0), 0.0, 1.0/sqrt(2.0))
  assert list(exp.beam.get_s0()) == list(s0_vec)
  assert list(exp.goniometer.get_rotation_axis()) == list((0.0, 0.0, 1.0))
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
  assert approx_equal(list(rotated_vectors), list(flex.vec3_double([(0.0, 0.0, 1.0),
    (0.0, 1.0, 0.0), (-1.0, 0.0, 0.0), (-1.0, 0.0, 1.0)])))

def test_sph_harm_table():
  """Simple test for the spherical harmonic table, constructing the table step
  by step, and verifying the values of a few easy-to-calculate entries.
  This also acts as a test for the calc_theta_phi function as well."""
  rt, exp = generate_input_vectors()
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
