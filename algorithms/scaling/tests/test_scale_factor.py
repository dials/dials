"""
Tests for scale components module.
"""
from math import exp
from dials.array_family import flex
from model.components.scale_components import SHScaleComponent, \
  SingleBScaleFactor, SingleScaleFactor, ScaleComponentBase
from model.components.smooth_scale_components import SmoothScaleComponent1D,\
  SmoothBScaleComponent1D, SmoothScaleComponent2D, SmoothScaleComponent3D
import pytest
from scitbx import sparse


def test_ScaleComponentBase():
  """Test for the ScaleComponentBase class."""

  class base_SF_filler(ScaleComponentBase):
    """Subclass to fill in the abstract method."""
    def update_reflection_data(self, reflection_table, selection=None):
      pass
    def calculate_scales_and_derivatives(self):
      pass
    def calculate_scales_derivatives_curvatures(self):
      pass

  # Test initialisation with no parameter esds.
  base_SF = base_SF_filler(flex.double([1.0] * 3))
  assert base_SF.n_params == 3
  assert base_SF.parameter_esds is None
  assert list(base_SF.parameters) == list(flex.double([1.0, 1.0, 1.0]))

  # Test updating of parameters
  base_SF.parameters = flex.double([2.0, 2.0, 2.0])
  assert list(base_SF.parameters) == list(flex.double([2.0, 2.0, 2.0]))
  with pytest.raises(AssertionError):
    # Try to change the number of parameters - should fail
    base_SF.parameters = flex.double([2.0, 2.0, 2.0, 2.0])

  # Test setting of inverse scales and updating to a list of different length.
  base_SF.inverse_scales = flex.double([1.0, 1.0])
  assert list(base_SF.inverse_scales) == list(flex.double([1.0, 1.0]))
  base_SF.inverse_scales = flex.double([1.0, 1.0, 1.0, 1.0])
  assert list(base_SF.inverse_scales) == list(flex.double([1.0, 1.0, 1.0, 1.0]))

  # Test setting of var_cov matrix.
  assert base_SF.var_cov_matrix is None
  base_SF.var_cov_matrix = [1.0, 0.0, 0.0]
  assert base_SF.var_cov_matrix == [1.0, 0.0, 0.0]

  assert hasattr(base_SF, '_derivatives')
  assert hasattr(base_SF, '_curvatures')

  base_SF = base_SF_filler(flex.double(3, 1.0), flex.double(3, 0.1))
  assert list(base_SF.parameter_esds) == list(flex.double(3, 0.1))

def test_SingleScaleFactor():
  """Test for SingleScaleFactor class."""
  KSF = SingleScaleFactor(flex.double([2.0]))
  assert KSF.n_params == 1
  assert list(KSF.parameters) == list(flex.double([2.0]))
  rt = flex.reflection_table()
  rt['d'] = flex.double([1.0, 1.0])
  KSF.update_reflection_data(rt)
  assert KSF.n_refl == 2
  KSF.calculate_scales_and_derivatives()
  assert list(KSF.inverse_scales) == list(flex.double([2.0, 2.0]))
  assert KSF.derivatives[0, 0] == 1
  assert KSF.derivatives[1, 0] == 1
  KSF.calculate_scales_derivatives_curvatures()
  assert KSF.curvatures.non_zeroes == 0
  KSF.update_reflection_data(rt, flex.bool([True, False])) # Test selection.
  assert KSF.n_refl == 1

def test_SingleBScaleFactor():
  """Test forSingleBScaleFactor class."""
  BSF = SingleBScaleFactor(flex.double([0.0]))
  assert BSF.n_params == 1
  assert list(BSF.parameters) == list(flex.double([0.0]))
  rt = flex.reflection_table()
  rt['d'] = flex.double([1.0, 1.0])
  BSF.update_reflection_data(rt)
  assert BSF.n_refl == 2
  assert list(BSF.d_values) == list(flex.double([1.0, 1.0]))
  BSF.calculate_scales_and_derivatives()
  assert list(BSF.inverse_scales) == list(flex.double([1.0, 1.0]))
  assert BSF.derivatives[0, 0] == 0.5
  assert BSF.derivatives[1, 0] == 0.5
  BSF.calculate_scales_derivatives_curvatures()
  assert BSF.curvatures[0, 0] == 0.25
  assert BSF.curvatures[1, 0] == 0.25
  BSF.update_reflection_data(rt, flex.bool([True, False])) # Test selection.
  assert BSF.n_refl == 1

def test_SHScalefactor():
  """Test the spherical harmonic absorption component."""
  initial_param = 0.1
  initial_val = 0.2

  SF = SHScaleComponent(flex.double([initial_param] * 3))
  assert SF.n_params == 3
  assert list(SF.parameters) == list(flex.double([initial_param]*3))

  # Test functionality just by setting sph_harm_table directly and calling
  # update_reflection_data to initialise the harmonic values.
  harmonic_values = sparse.matrix(1, 3)
  harmonic_values[0, 0] = initial_val
  harmonic_values[0, 1] = initial_val
  harmonic_values[0, 2] = initial_val
  SF.sph_harm_table = harmonic_values
  SF.update_reflection_data(None, flex.bool([True]))
  assert SF.harmonic_values[0, 0] == initial_val
  assert SF.harmonic_values[0, 1] == initial_val
  assert SF.harmonic_values[0, 2] == initial_val
  assert list(SF.inverse_scales) == [1.0 + (3.0 * initial_val * initial_param)]
  assert SF.derivatives[0, 0] == initial_val
  assert SF.derivatives[0, 1] == initial_val
  assert SF.derivatives[0, 2] == initial_val
  SF.calculate_scales_derivatives_curvatures()
  assert SF.curvatures.non_zeroes == 0

def test_SmoothScaleFactor1D():
  """Test for the gaussian smoothed 1D scalefactor class."""
  SF = SmoothScaleComponent1D(flex.double(5, 1.1), col_name='norm_rot')
  assert SF.n_params == 5
  assert list(SF.parameters) == list(flex.double([1.1, 1.1, 1.1, 1.1, 1.1]))
  rt = flex.reflection_table()
  rt['norm_rot'] = flex.double([0.5, 1.0, 2.5, 0.0])
  SF.update_reflection_data(rt)
  assert list(SF.normalised_values) == list(flex.double([0.5, 1.0, 2.5, 0.0]))
  assert list(SF.inverse_scales) == list(flex.double([1.0, 1.0, 1.0, 1.0]))
  SF._smoother.set_smoothing(4, 1.0)
  SF.calculate_scales()
  assert (list(abs(SF.inverse_scales - flex.double([1.1, 1.1, 1.1, 1.1]))) <
    list(flex.double([1e-7]*len(SF.inverse_scales))))
  SF.calculate_scales_and_derivatives()
  assert abs((SF.derivatives[0, 0]/SF.derivatives[0, 1]) - exp(-1.0)/exp(0.0)) < 1e-6
  T = SF.derivatives.transpose()
  assert sum(list(T[:, 0].as_dense_vector())) == 1.0 #should always be 1.0
  assert sum(list(T[:, 1].as_dense_vector())) == 1.0
  assert sum(list(T[:, 2].as_dense_vector())) == 1.0
  assert SF.derivatives[1, 1] == SF.derivatives[1, 2]
  assert SF.derivatives[1, 0] == SF.derivatives[1, 3]

def test_SmoothBScaleFactor1D():
  'test for a gaussian smoothed 1D scalefactor object'
  SF = SmoothBScaleComponent1D(flex.double(5, 0.0), col_name='norm_rot')
  assert SF.n_params == 5
  assert list(SF.parameters) == list(flex.double(5, 0.0))
  rt = flex.reflection_table()
  rt['norm_rot'] = flex.double([0.5, 1.0, 2.5, 0.0])
  rt['d'] = flex.double([1.0, 1.0, 1.0, 1.0])
  SF.update_reflection_data(rt)
  assert list(SF.normalised_values) == list(flex.double([0.5, 1.0, 2.5, 0.0]))
  assert list(SF.d_values) == list(flex.double([1.0, 1.0, 1.0, 1.0]))
  assert list(SF.inverse_scales) == list(flex.double([1.0, 1.0, 1.0, 1.0]))
  SF._smoother.set_smoothing(4, 1.0)
  SF.calculate_scales()
  assert (list(abs(SF.inverse_scales - flex.double([1.0, 1.0, 1.0, 1.0]))) <
    list(flex.double([1e-7]*len(SF.inverse_scales))))
  SF.calculate_scales_and_derivatives()
  assert abs((SF.derivatives[0, 0]/SF.derivatives[0, 1]) - exp(-1.0)/exp(0.0)) < 1e-6
  T = SF.derivatives.transpose()
  assert sum(list(T[:, 0].as_dense_vector())) == 0.5 #value depends on d
  assert sum(list(T[:, 1].as_dense_vector())) == 0.5
  assert sum(list(T[:, 2].as_dense_vector())) == 0.5
  assert SF.derivatives[1, 1] == SF.derivatives[1, 2]
  assert SF.derivatives[1, 0] == SF.derivatives[1, 3]
