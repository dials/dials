"""
Tests for scale components module.
"""
from dials.array_family import flex
from model.components.scale_components import SHScaleComponent, \
  SingleBScaleFactor, SingleScaleFactor, ScaleComponentBase
from model.components.smooth_scale_components import SmoothScaleComponent1D,\
  SmoothBScaleComponent1D, SmoothScaleComponent2D, SmoothScaleComponent3D
import pytest

def test_base_scalefactor():
  'test for abstract scalefactor object'

  class base_SF_filler(ScaleComponentBase):
    'subclass ScaleFactor to fill in the abstract method'
    def update_reflection_data(self, reflection_table, selection=None):
      pass
    def calculate_scales_and_derivatives(self):
      pass

  base_SF = base_SF_filler(flex.double([1.0] * 3))
  assert base_SF.n_params == 3
  assert list(base_SF.parameters) == list(flex.double([1.0, 1.0, 1.0]))

  base_SF.parameters = flex.double([2.0, 2.0, 2.0])
  assert list(base_SF.parameters) == list(flex.double([2.0, 2.0, 2.0]))
  with pytest.raises(AssertionError):
    #try to change the number of parameters - should fail
    base_SF.parameters = flex.double([2.0, 2.0, 2.0, 2.0])

  base_SF.inverse_scales = flex.double([1.0, 1.0])
  assert list(base_SF.inverse_scales) == list(flex.double([1.0, 1.0]))
  base_SF.inverse_scales = flex.double([1.0, 1.0, 1.0, 1.0])
  assert list(base_SF.inverse_scales) == list(flex.double([1.0, 1.0, 1.0, 1.0]))

def test_KScaleFactor():
  'test for simple global K scalefactor object'
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

def test_BScaleFactor():
  'test for simple global B scalefactor object'
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

def test_SmoothScaleFactor1D():
  'test for a gaussian smoothed 1D scalefactor object'
  from math import exp
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
  from math import exp
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


'''def test_SHScalefactor():
  """Test the spherical harmonic absorption component."""
  initial_param = 0.1
  initial_val = 0.2

  rt = flex.reflection_table()
  from math import sqrt, pi
  s1_vec = (1.0/sqrt(2.0), 0.0, 1.0/sqrt(2.0))
  s0_vec = (1.0, 0.0, 0.0)
  rt['s1'] = flex.vec3_double([s1_vec, s1_vec, s1_vec])
  rt['phi'] = flex.double([0.0, 45.0, 90.0])
  exp = generated_single_exp()[0]
  param = generated_param()

  SF = SHScaleComponent(flex.double([initial_param] * 3))
  assert SF.n_params == 3
  assert list(SF.parameters) == list(flex.double([initial_param]*3))

  SF.configure_reflection_table(rt, exp, param)

  harmonic_values = sparse.matrix(1, 3)
  harmonic_values[0, 0] = initial_val
  harmonic_values[0, 1] = initial_val
  harmonic_values[0, 2] = initial_val
  SF.update_reflection_data(harmonic_values)
  assert SF.harmonic_values[0, 0] == initial_val
  assert SF.harmonic_values[0, 1] == initial_val
  assert SF.harmonic_values[0, 2] == initial_val
  assert list(SF.inverse_scales) == [1.0 + (3.0 * initial_val * initial_param)]
  assert SF.derivatives[0, 0] == initial_val
  assert SF.derivatives[0, 1] == initial_val
  assert SF.derivatives[0, 2] == initial_val'''
