"""
Test the output of targeted scaling - by calling the scale_against_target
scaling_library function and by directly invoking the perform method scaling
of a TargetScaler.
"""

from math import log
import pytest
from dials.array_family import flex
from dials.util.options import OptionParser
from libtbx import phil
from libtbx.test_utils import approx_equal
from dxtbx.model import Crystal, Experiment, ExperimentList
from dials.algorithms.scaling.scaler_factory import TargetScalerFactory
from dials.algorithms.scaling.scaling_library import scale_against_target

def generated_target_refl():
  """Test target reflection table."""
  reflections = flex.reflection_table()
  reflections['intensity.prf.value'] = flex.double([1.0, 10.0, 1.0])
  reflections['intensity.prf.variance'] = flex.double([1.0, 10.0, 1.0])
  reflections['miller_index'] = flex.miller_index([(1, 0, 0), (0, 0, 1),
    (1, 0, 0)]) #don't change
  reflections['d'] = flex.double([1.0, 2.0, 1.0])
  reflections['partiality'] = flex.double([1.0, 1.0, 1.0])
  reflections.set_flags(flex.bool([True, True, True]),
    reflections.flags.integrated)
  return reflections

def generated_refl_to_scale():
  """Generate a reflection table for targeted scaling."""
  reflections = flex.reflection_table()
  reflections['intensity.prf.value'] = flex.double([2.0, 5.0, 2.0, 1.0])
  reflections['intensity.prf.variance'] = flex.double([2.0, 5.0, 2.0, 1.0])
  reflections['miller_index'] = flex.miller_index([(1, 0, 0), (0, 0, 1),
    (1, 0, 0), (10, 0, 0)]) #don't change
  reflections['d'] = flex.double([1.0, 2.0, 1.0, (4.0/3.0)**0.5])
  reflections['partiality'] = flex.double([1.0, 1.0, 1.0, 1.0])
  reflections.set_flags(flex.bool([True, True, True, True]),
    reflections.flags.integrated)
  return reflections

@pytest.fixture
def test_target_refl():
  """Return the target reflection table."""
  return generated_target_refl()

@pytest.fixture()
def test_refl_to_scale():
  """Return the reflection table to scale."""
  return generated_refl_to_scale()

@pytest.fixture
def test_exp():
  """Test experiments object."""
  experiments = ExperimentList()
  exp_dict = {"__id__" : "crystal", "real_space_a": [1.0, 0.0, 0.0],
              "real_space_b": [0.0, 1.0, 0.0], "real_space_c": [0.0, 0.0, 2.0],
              "space_group_hall_symbol": " C 2y"}
  crystal = Crystal.from_dict(exp_dict)
  experiments.append(Experiment(crystal=crystal))
  return experiments

@pytest.fixture
def KB_test_param():
  """Generate a params phil scope."""
  phil_scope = phil.parse('''
      include scope dials.algorithms.scaling.scaling_options.phil_scope
      include scope dials.algorithms.scaling.scaling_refiner.scaling_refinery_phil_str
  ''', process_includes=True)

  optionparser = OptionParser(phil=phil_scope, check_format=False)
  parameters, _ = optionparser.parse_args(args=[], quick_parse=True,
    show_diff_phil=False)
  parameters.__inject__('model', 'KB')
  return parameters


def test_scale_against_target(test_refl_to_scale, test_target_refl,
  test_exp, KB_test_param):
  """Integration testing of the scale_against_target library function/targeted
  scaling."""
  # Based on input - have Target Ih of 1.0, 10.0 and d of 1.0, 2.0. Then refl to
  # target have I's of 2 and 5, (with the same ds). Therefore for a KB model the
  # problem can be minimised exactly by solving the equations:
  # 2 = K * exp(B/2)
  # 1/2 = K * exp(B/8)
  # Solving these gives the form tested for at the end of this test.
  target_reflections = test_target_refl
  reflections = test_refl_to_scale
  target_experiments = test_exp
  experiments = test_exp
  scaled_reflections = scale_against_target(reflections, experiments,
    target_reflections, target_experiments)
  assert approx_equal(list(scaled_reflections['inverse_scale_factor']),
    [2.0, 0.5, 2.0, 2.0 * (4.0 **(-1.0/3.0))])

  # Repeat the test but calling the TargetScaler directly, to allow inspection
  # of the model components.
  targetscaler = TargetScalerFactory.create(KB_test_param, experiments,
    [target_reflections, reflections], [True, False])
  targetscaler.perform_scaling()
  assert approx_equal(list(targetscaler.unscaled_scalers[0].components[
    'scale'].parameters), [(4.0**(-1.0/3.0))/2.0])
  assert approx_equal(list(targetscaler.unscaled_scalers[0].components[
    'decay'].parameters), [(log(4.0) * 8.0/3.0)])
