"""
Tests for the error model.
"""
from math import sqrt
import pytest
from libtbx.test_utils import approx_equal
from libtbx.utils import Sorry
from dials.algorithms.scaling.error_model.error_model import \
  get_error_model, BasicErrorModel
from dials.algorithms.scaling.error_model.error_model_target import ErrorModelTarget
from dials.algorithms.scaling.Ih_table import IhTable
from dials.array_family import flex
from cctbx.sgtbx import space_group

@pytest.fixture()
def large_reflection_table():
  """Create a reflection table."""
  return generate_refl_1()

@pytest.fixture(scope='module')
def test_sg():
  """Create a space group object."""
  return space_group("P 1")

def generate_refl_1():
  """Generate a test reflection table. Note tha the variance values are chosen
  as the 'True' Ih_values, which would be found if unity weights were chosen
  in this example."""
  reflections = flex.reflection_table()
  reflections['intensity'] = flex.double([1.0, 2.0, 3.0, 4.0, 5.0,
    6.0, 7.0, 8.0, 9.0, 10.0])
  reflections['inverse_scale_factor'] = flex.double(10, 1.0)
  reflections['variance'] = flex.double([1.0, 2.0, 3.0, 4.0, 5.0,
    6.0, 7.0, 8.0, 9.0, 10.0])
  reflections['miller_index'] = flex.miller_index([(1, 0, 0), (0, 0, 1),
    (-3, 0, 0), (0, 2, 0), (1, 0, 0), (0, 0, -20), (0, 0, 2), (10, 0, 0),
    (0, 0, 10), (1, 0, 0)])
  reflections.set_flags(flex.bool(10, True), reflections.flags.integrated)
  return reflections

def test_errormodel(large_reflection_table, test_sg):
  """Test the initialisation and methods of the error model."""

  # first test get_error_model helper function.
  with pytest.raises(Sorry):
    em = get_error_model('bad')
  em = get_error_model('basic')
  Ih_table = IhTable([large_reflection_table], test_sg, nblocks=1)
  block = Ih_table.blocked_data_list[0]
  error_model = em(block, n_bins=10)
  assert error_model.summation_matrix[0, 1] == 1
  assert error_model.summation_matrix[1, 6] == 1
  assert error_model.summation_matrix[2, 8] == 1
  assert error_model.summation_matrix[3, 5] == 1
  assert error_model.summation_matrix[4, 3] == 1
  assert error_model.summation_matrix[5, 0] == 1
  assert error_model.summation_matrix[6, 4] == 1
  assert error_model.summation_matrix[7, 9] == 1
  assert error_model.summation_matrix[8, 2] == 1
  assert error_model.summation_matrix[9, 7] == 1
  assert error_model.summation_matrix.non_zeroes == large_reflection_table.size()
  assert error_model.bin_counts == flex.double(large_reflection_table.size(), 1)
  assert list(error_model.n_h) == list(block.calc_nh())

  # Test calc sigmaprime
  x0 = 1.0
  x1 = 0.1
  error_model.sigmaprime = error_model.calc_sigmaprime([x0, x1])
  assert error_model.sigmaprime == x0 * ((block.variances +
    ((x1 * block.intensities)**2))**0.5)

  # Test calc delta_hl
  error_model.sigmaprime = error_model.calc_sigmaprime([1.0, 0.0]) #Reset
  # Calculate example for three elements, with intensities 1, 5 and 10 and
  # variances 1, 5 and 10 using he formula
  # delta_hl = sqrt(n_h - 1 / n_h) * (Ihl/ghl - Ih) / sigmaprime
  error_model.delta_hl = error_model.calc_deltahl()
  expected_deltas = [0.0, 0.0, 0.0, 0.0, 0.0, (-17.0/13.0) * sqrt(2.0/3.0),
    (7.0/13.0) * sqrt(10.0/3.0), (10.0/13.0) * sqrt(20.0/3.0), 0.0, 0.0]
  assert approx_equal(list(error_model.delta_hl), expected_deltas)

  # Test bin variance calculation on example with fewer bins.
  error_model = BasicErrorModel(block, n_bins=2)

  assert error_model.summation_matrix[0, 0] == 1
  assert error_model.summation_matrix[1, 1] == 1
  assert error_model.summation_matrix[2, 1] == 1
  assert error_model.summation_matrix[3, 1] == 1
  assert error_model.summation_matrix[4, 0] == 1
  assert error_model.summation_matrix[5, 0] == 1
  assert error_model.summation_matrix[6, 0] == 1
  assert error_model.summation_matrix[7, 1] == 1
  assert error_model.summation_matrix[8, 0] == 1
  assert error_model.summation_matrix[9, 1] == 1

  error_model.sigmaprime = error_model.calc_sigmaprime([1.0, 0.0])
  error_model.delta_hl = error_model.calc_deltahl()
  bin_vars = error_model.calculate_bin_variances()
  mu_0 = sqrt(2.0/3.0) * ((7.0 * sqrt(5.0)) - 17.0)/ 65.0
  a = ((3.0 * mu_0**2) + (((7.0/13.0) * sqrt(10.0/3.0)) - mu_0)**2 +
    (((-17.0/13.0) * sqrt(2.0/3.0)) - mu_0)**2)/5.0
  expected_bin_vars = [a, 320.0/507.0]
  assert approx_equal(list(bin_vars), expected_bin_vars)

  # Test updating call
  error_model = BasicErrorModel(block, n_bins=2)
  error_model.update_for_minimisation([1.0, 0.0])
  assert approx_equal(list(error_model.bin_variances), expected_bin_vars)

def test_error_model_target(large_reflection_table, test_sg):
  """Test the error model target."""
  Ih_table = IhTable([large_reflection_table], test_sg, nblocks=1)
  block = Ih_table.blocked_data_list[0]
  error_model = BasicErrorModel(block, n_bins=2)
  error_model.update_for_minimisation([1.0, 0.05])
  target = ErrorModelTarget(error_model)
  # Test residual calculation
  residuals = target.calculate_residuals()
  assert residuals == (flex.double(2, 1.0) - error_model.bin_variances)**2

  # Test gradient calculation against finite differences.
  gradients = target.calculate_gradients()
  gradient_fd = calculate_gradient_fd(target)
  assert approx_equal(gradients, gradient_fd)

  # Test the method calls
  r, g = target.compute_functional_gradients()
  assert r == residuals
  assert list(gradients) == list(g)
  r, g, c = target.compute_functional_gradients_and_curvatures()
  assert r == residuals
  assert list(gradients) == list(g)
  assert c is None

def calculate_gradient_fd(target):
  """Calculate gradient array with finite difference approach."""
  delta = 1.0e-6
  gradients = flex.double([0.0] * len(target.x))
  #iterate over parameters, varying one at a time and calculating the gradient
  for i in range(len(target.x)):
    target.x[i] -= 0.5 * delta
    target.predict()
    R_low = target.calculate_residuals()
    target.x[i] += delta
    target.predict()
    R_upper = target.calculate_residuals()
    target.x[i] -= 0.5 * delta
    target.predict()
    gradients[i] = (flex.sum(R_upper) - flex.sum(R_low)) / delta
  return gradients
