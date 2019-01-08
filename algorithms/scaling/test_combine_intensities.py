from math import sqrt, pi
import pytest
import numpy as np
from mock import Mock
from dxtbx.model import Experiment, Crystal
from libtbx.test_utils import approx_equal
from dials.array_family import flex
from dials.algorithms.scaling.combine_intensities import combine_intensities, optimise_intensity_combination

@pytest.fixture(scope='module')
def test_exp_P1():
  """Create a mock experiments object."""
  exp = Experiment()
  exp_dict = {"__id__" : "crystal", "real_space_a": [1.0, 0.0, 0.0],
              "real_space_b": [0.0, 1.0, 0.0], "real_space_c": [0.0, 0.0, 1.0],
              "space_group_hall_symbol": " P 1"}
  crystal = Crystal.from_dict(exp_dict)
  exp.crystal = crystal
  return exp

def generate_simple_table(prf=True):
  """Generate a reflection table for testing intensity combination.
  The numbers are contrived to make sum intensities agree well at high
  intensity but terribly at low and vice versa for profile intensities."""
  reflections = flex.reflection_table()
  reflections['miller_index'] = flex.miller_index([
    (0, 0, 1), (0, 0, 1), (0, 0, 1), (0, 0, 1), (0, 0, 1),
    (0, 0, 2), (0, 0, 2), (0, 0, 2), (0, 0, 2), (0, 0, 2),
    (0, 0, 3), (0, 0, 3), (0, 0, 3), (0, 0, 3), (0, 0, 3),
    (0, 0, 4), (0, 0, 4), (0, 0, 4), (0, 0, 4), (0, 0, 4),
    (0, 0, 5), (0, 0, 5), (0, 0, 5), (0, 0, 5), (0, 0, 5)])
  reflections['inverse_scale_factor'] = flex.double(25, 1.0)
  #Contrive an example that should give the best cc12 when combined.
  #make sum intensities agree well at high intensity but terribly at low
  # and vice versa for profile intensities.
  #profile less consistent at high intensity here

  #sumless consistent at low intensity here
  reflections['intensity.sum.value'] = flex.double([
    10000.0, 11000.0, 9000.0, 8000.0, 12000.0,
    500.0, 5600.0, 5500.0, 2000.0, 6000.0,
    100.0, 50.0, 150.0, 75.0, 125.0,
    30.0, 10.0, 2.0, 35.0, 79.0,
    1.0, 10.0, 20.0, 10.0, 5.0])
  reflections['intensity.sum.variance'] = flex.double(
    [10000]*5 + [5000]*5 + [100]*5 + [30]*5 + [10]*5)
  reflections.set_flags(flex.bool(25, False), reflections.flags.outlier_in_scaling)
  reflections.set_flags(flex.bool(25, True), reflections.flags.integrated)
  if prf:
    reflections['intensity.prf.value'] = flex.double([
      10000.0, 16000.0, 12000.0, 6000.0, 9000.0,
      5000.0, 2000.0, 1500.0, 1300.0, 9000.0,
      100.0, 80.0, 120.0, 90.0, 100.0,
      30.0, 40.0, 50.0, 30.0, 30.0,
      10.0, 12.0, 9.0, 8.0, 10.0])
    reflections['intensity.prf.variance'] = flex.double(
      [10000]*5 + [5000]*5 + [100]*5 + [30]*5 + [10]*5)
  return reflections

def test_optimise_intensity_combination(test_exp_P1):
  """Test optimise_intensity_combination function."""
  r1 = flex.reflection_table()
  r2 = flex.reflection_table()

  r1['intensity.prf.value'] = flex.double(range(1, 101))
  r1['intensity.sum.value'] = flex.double(range(2, 102))
  r1['intensity.prf.variance'] = flex.double(range(1, 101))
  r1['intensity.sum.variance'] = flex.double(range(2, 102))
  r1['inverse_scale_factor'] = flex.double(100, 1)
  r1['miller_index'] = flex.miller_index([(1,0,0)]*100)

  r2['intensity.prf.value'] = flex.double(range(1, 101))
  r2['intensity.sum.value'] = flex.double(range(2, 102))
  r2['intensity.prf.variance'] = flex.double(range(1, 101))
  r2['intensity.sum.variance'] = flex.double(range(2, 102))
  r2['inverse_scale_factor'] = flex.double(100, 1)
  r2['miller_index'] = flex.miller_index([(1,0,0)]*100)

  r1.set_flags(flex.bool(100, True), r1.flags.integrated_sum)
  r1.set_flags(flex.bool([True]*98 + [False]*2), r1.flags.integrated_prf)
  r2.set_flags(flex.bool(100, True), r2.flags.integrated_sum)
  r2.set_flags(flex.bool([True]*98 + [False]*2), r2.flags.integrated_prf)

  _ = optimise_intensity_combination([r1], test_exp_P1)

  _ = optimise_intensity_combination([r1, r2], test_exp_P1)

def test_combine_intensities(test_exp_P1):
  """Test the combine intensities function for a single dataset"""
  reflections = generate_simple_table()
  Imid = optimise_intensity_combination([reflections], test_exp_P1)
  reflections = combine_intensities(reflections, Imid)
  #reflections_list, results = combine_intensities([reflections], test_exp_P1)
  # Imid being 1200.0 should be best for this contrived example
  assert Imid == 1200.0

  #Due to nature of crossover, just require 2% tolerance for this example
  assert list(reflections['intensity'][0:5]) == pytest.approx(list(
    reflections['intensity.sum.value'][0:5]), rel=2e-2)
  assert list(reflections['intensity'][20:25]) == pytest.approx(list(
    reflections['intensity.prf.value'][20:25]), rel=2e-2)

def test_combine_intensities_multi_dataset(test_exp_P1):
  """Test the combine intensities function for multiple datasets"""
  r1 = generate_simple_table()
  r1['partiality'] = flex.double(25, 1.0)
  r2 = generate_simple_table(prf=False)
  Imid = optimise_intensity_combination([r1, r2], test_exp_P1)
  assert pytest.approx(Imid) == 1200.0

  r1 = generate_simple_table()
  r1['partiality'] = flex.double(25, 1.0)
  r1 = combine_intensities(r1, 0)
  assert list(r1['intensity']) == list(r1['intensity.prf.value'])
  r1 = combine_intensities(r1, 1)
  assert list(r1['intensity']) == list(r1['intensity.sum.value'])

  r2 = generate_simple_table(prf=False)
  r2 = combine_intensities(r2, 1)
  assert list(r2['intensity']) == list(r2['intensity.sum.value'])

  r1 = generate_simple_table(prf=False)
  r2 = generate_simple_table(prf=False)
  Imid = optimise_intensity_combination([r1, r2], test_exp_P1)
  r1 = combine_intensities(r1, Imid)
  assert list(r1['intensity']) == list(r1['intensity.sum.value'])
