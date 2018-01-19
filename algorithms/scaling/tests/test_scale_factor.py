'''
tests for scalefactor objects
'''
from dials.array_family import flex
from scale_factor import (ScaleFactor, KScaleFactor, BScaleFactor,
  SmoothScaleFactor, SmoothScaleFactor1D, SmoothBScaleFactor1D, SHScaleFactor)
import pytest
from scitbx import sparse

def test_base_scalefactor():
  'test for abstract scalefactor object'

  class base_SF_filler(ScaleFactor):
    'subclass ScaleFactor to fill in the abstract method'
    def update_reflection_data(self):
      pass

  base_SF = base_SF_filler(1.0, 3)
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
  KSF = KScaleFactor(2.0)
  assert KSF.n_params == 1
  assert list(KSF.parameters) == list(flex.double([2.0]))
  KSF.update_reflection_data(n_refl=2)
  assert KSF.n_refl == 2
  KSF.calculate_scales_and_derivatives()
  assert list(KSF.inverse_scales) == list(flex.double([2.0, 2.0]))
  assert KSF.derivatives[0, 0] == 1
  assert KSF.derivatives[1, 0] == 1

def test_BScaleFactor():
  'test for simple global B scalefactor object'
  BSF = BScaleFactor(0.0)
  assert BSF.n_params == 1
  assert list(BSF.parameters) == list(flex.double([0.0]))
  BSF.update_reflection_data(dvalues=flex.double([1.0, 1.0]))
  assert BSF.n_refl == 2
  assert list(BSF.d_values) == list(flex.double([1.0, 1.0]))
  BSF.calculate_scales_and_derivatives()
  assert list(BSF.inverse_scales) == list(flex.double([1.0, 1.0]))
  assert BSF.derivatives[0, 0] == 0.5
  assert BSF.derivatives[1, 0] == 0.5

def test_SmoothScaleFactor1D():
  'test for a gaussian smoothed 1D scalefactor object'
  from math import exp
  SF = SmoothScaleFactor1D(1.1, 5)
  assert SF.n_params == 5
  assert list(SF.parameters) == list(flex.double([1.1, 1.1, 1.1, 1.1, 1.1]))
  SF.update_reflection_data(normalised_values=[0.5, 1.0, 2.5])
  assert list(SF.normalised_values) == list(flex.double([0.5, 1.0, 2.5]))
  assert list(SF.inverse_scales) == list(flex.double([1.0, 1.0, 1.0]))
  SF._smoother.set_smoothing(4, 1.0)
  SF.calculate_scales()
  assert (list(abs(SF.inverse_scales - flex.double([1.1, 1.1, 1.1]))) <
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
  SF = SmoothBScaleFactor1D(0.0, 5)
  assert SF.n_params == 5
  assert list(SF.parameters) == list(flex.double([0.0]*5))
  SF.update_reflection_data(normalised_values=[0.5, 1.0, 2.5],
    dvalues=flex.double([1.0, 1.0, 1.0]))
  assert list(SF.normalised_values) == list(flex.double([0.5, 1.0, 2.5]))
  assert list(SF.d_values) == list(flex.double([1.0, 1.0, 1.0]))
  assert list(SF.inverse_scales) == list(flex.double([1.0, 1.0, 1.0]))
  SF._smoother.set_smoothing(4, 1.0)
  SF.calculate_scales()
  assert (list(abs(SF.inverse_scales - flex.double([1.0, 1.0, 1.0]))) <
    list(flex.double([1e-7]*len(SF.inverse_scales))))
  SF.calculate_scales_and_derivatives()
  assert abs((SF.derivatives[0, 0]/SF.derivatives[0, 1]) - exp(-1.0)/exp(0.0)) < 1e-6
  T = SF.derivatives.transpose()
  assert sum(list(T[:, 0].as_dense_vector())) == 0.5 #value depends on d
  assert sum(list(T[:, 1].as_dense_vector())) == 0.5
  assert sum(list(T[:, 2].as_dense_vector())) == 0.5
  assert SF.derivatives[1, 1] == SF.derivatives[1, 2]
  assert SF.derivatives[1, 0] == SF.derivatives[1, 3]

def test_SHScalefactor():
  initial_param = 0.1
  initial_val = 0.2
  SF = SHScaleFactor(initial_param, 3)
  assert SF.n_params == 3
  assert list(SF.parameters) == list(flex.double([initial_param]*3))
  harmonic_values = sparse.matrix(1, 3)
  harmonic_values[0, 0] = initial_val
  harmonic_values[0, 1] = initial_val
  harmonic_values[0, 2] = initial_val
  SF.update_reflection_data(harmonic_values)
  assert SF.harmonic_values[0, 0] == initial_val
  assert SF.harmonic_values[0, 1] == initial_val
  assert SF.harmonic_values[0, 2] == initial_val
  assert list(SF.inverse_scales) == [1.0 + (3.0 * initial_val * initial_param)] #is 1.0 + 3 x param x harm_val
  assert SF.derivatives[0, 0] == initial_val
  assert SF.derivatives[0, 1] == initial_val
  assert SF.derivatives[0, 2] == initial_val
