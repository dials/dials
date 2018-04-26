from math import exp
from libtbx.test_utils import approx_equal
from dials.array_family import flex
from dials.algorithms.scaling.model.components.smooth_scale_components import \
  GaussianSmoother2D, GaussianSmoother3D

def test_2DGaussianSmoother():
  """Test the behvaiour of the 2D Gaussian Smoother"""
  GS2D = GaussianSmoother2D([0, 1], 1, [0, 4], 4)
  assert GS2D.num_x_values() == 2
  assert GS2D.num_y_values() == 6
  assert GS2D.num_x_average() == 2
  assert GS2D.num_y_average() == 3
  parameters = flex.double([1.0, 1.0, 2.0, 2.0, 3.0, 3.0, 4.0, 4.0, 5.0, 5.0,
    6.0, 6.0])
  x = flex.double([0.5, 0.25])
  y = flex.double([1.0, 0.75])
  value, weight, sumw = GS2D.multi_value_weight(x, y, parameters)
  assert weight.non_zeroes == 6 * 2
  assert GS2D.sigma() == 0.65
  # Calculate tand verify the expected value at the first position.
  expected_sumw = ((2.0*exp(-2.5/(0.65**2)))+ 4.0*exp(-0.5/(0.65**2)))
  assert approx_equal(expected_sumw, sumw[0])
  expected_value_numerator = ((10.0 * exp(-0.5/(0.65**2))) +
    (2.0 * exp(-2.5/(0.65**2))))
  assert approx_equal(value[0], expected_value_numerator / expected_sumw)

  # Do the same calcualtion but calling the single value_weight function.
  x = 0.5
  y = 1.0
  value, weight, sumw = GS2D.value_weight(x, y, parameters)
  assert weight.non_zeroes == 6
  assert approx_equal(value, expected_value_numerator / expected_sumw)

  #Test again for another small number of params
  GS2D = GaussianSmoother2D([0, 2], 2, [0, 4], 4)
  assert GS2D.num_x_values() == 3
  assert GS2D.num_y_values() == 6
  assert GS2D.num_x_average() == 3
  assert GS2D.num_y_average() == 3

  #Test again for another number of params
  GS2D = GaussianSmoother2D([0, 4], 4, [0, 5], 5)
  assert GS2D.num_x_values() == 6
  assert GS2D.num_y_values() == 7
  assert GS2D.num_x_average() == 3
  assert GS2D.num_y_average() == 3

def test_3DGaussianSmoother():
  """Test the behvaiour of the 2D Gaussian Smoother"""
  GS3D = GaussianSmoother3D([0, 1], 1, [0, 3], 3, [0, 1], 1)
  assert GS3D.num_x_values() == 2
  assert GS3D.num_y_values() == 5
  assert GS3D.num_z_values() == 2
  assert GS3D.num_x_average() == 2
  assert GS3D.num_y_average() == 3
  assert GS3D.num_z_average() == 2
  parameters = flex.double([1.0, 1.0, 2.0, 2.0, 3.0, 3.0, 4.0, 4.0, 5.0, 5.0] * 2)
  x = flex.double([0.5, 0.25])
  y = flex.double([1.0, 0.75])
  z = flex.double([0.5, 0.25])
  value, weight, sumw = GS3D.multi_value_weight(x, y, z, parameters)
  assert weight.non_zeroes == 6 * 2 * 2
  assert GS3D.sigma() == 0.65

  # Calculate tand verify the expected value at the first position.
  expected_sumw = (4.0 * exp(-2.75/(0.65**2)))+ (8.0 * exp(-0.75/(0.65**2)))
  assert approx_equal(expected_sumw, sumw[0])
  expected_value_numerator = ((20.0 * exp(-0.75/(0.65**2))) +
    (4.0 * exp(-2.75/(0.65**2))))
  assert approx_equal(value[0], expected_value_numerator / expected_sumw)

  # Do the same calcualtion but calling the single value_weight function.
  x = 0.5
  y = 1.0
  z = 0.5
  value, weight, sumw = GS3D.value_weight(x, y, z, parameters)
  assert weight.non_zeroes == 6 * 2
  assert approx_equal(value, expected_value_numerator / expected_sumw)

  #Test again for another number of params
  GS3D = GaussianSmoother3D([0, 1], 1, [0, 2], 2, [0, 3], 3)
  assert GS3D.num_x_values() == 2
  assert GS3D.num_y_values() == 3
  assert GS3D.num_z_values() == 5
  assert GS3D.num_x_average() == 2
  assert GS3D.num_y_average() == 3
  assert GS3D.num_y_average() == 3