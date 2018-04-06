"""
Test for the basis function module.
"""
import numpy as np
import pytest
from dials.array_family import flex
from dials.util.options import OptionParser
from parameter_handler import scaling_active_parameter_manager
from libtbx import phil
#from libtbx.test_utils import approx_equal
from dxtbx.model.experiment_list import ExperimentList
from dxtbx.model import Crystal, Scan, Beam, Goniometer, Detector, Experiment
from dials.algorithms.scaling.model.scaling_model_factory import \
  create_scaling_model
from dials.algorithms.scaling.scaler_factory import create_scaler
from basis_functions import basis_function


def generated_refl():
  """Generate reflection table to test the basis function."""
  #these miller_idx/d_values don't make physical sense, but I didn't want to
  #have to write the tests for lots of reflections.
  reflections = flex.reflection_table()
  reflections['intensity.prf.value'] = flex.double([1.0, 10.0, 100.0])
  reflections['intensity.prf.variance'] = flex.double([1.0, 10.0, 100.0])
  reflections['miller_index'] = flex.miller_index([(1, 0, 0), (0, 0, 1),
    (1, 0, 0)]) #don't change
  reflections['d'] = flex.double([0.8, 2.0, 2.0]) #don't change
  reflections['lp'] = flex.double([1.0, 1.0, 1.0])
  reflections['dqe'] = flex.double([1.0, 1.0, 1.0])
  reflections['partiality'] = flex.double([1.0, 1.0, 1.0])
  reflections['xyzobs.px.value'] = flex.vec3_double([(0.0, 0.0, 0.0),
    (0.0, 0.0, 5.0), (0.0, 0.0, 10.0)])
  reflections['s1'] = flex.vec3_double([(0.0, 0.1, 1.0), (0.0, 0.1, 1.0),
    (0.0, 0.1, 1.0)])
  reflections.set_flags(flex.bool([True, True, True]),
    reflections.flags.integrated)
  return [reflections]

#@pytest.fixture(scope='module')
def generated_single_exp():
  """Generate an experiment object."""
  experiments = ExperimentList()
  exp_dict = {"__id__" : "crystal", "real_space_a": [1.0, 0.0, 0.0],
              "real_space_b": [0.0, 1.0, 0.0], "real_space_c": [0.0, 0.0, 2.0],
              "space_group_hall_symbol": " C 2y"}
  crystal = Crystal.from_dict(exp_dict)
  scan = Scan(image_range=[0, 90], oscillation=[0.0, 1.0])
  beam = Beam(s0=(0.0, 0.0, 1.01))
  goniometer = Goniometer((1.0, 0.0, 0.0))
  detector = Detector()
  experiments.append(Experiment(beam=beam, scan=scan, goniometer=goniometer,
    detector=detector, crystal=crystal))
  return experiments

@pytest.fixture(scope='module')
def generated_KB_param():
  """Generate the scaling phil param scope."""
  phil_scope = phil.parse('''
      include scope dials.algorithms.scaling.scaling_options.phil_scope
  ''', process_includes=True)

  optionparser = OptionParser(phil=phil_scope, check_format=False)
  parameters, _ = optionparser.parse_args(args=None, quick_parse=True,
    show_diff_phil=False)
  parameters.__inject__('model', 'KB')
  return parameters



def test_basis_function(generated_KB_param):
  """Test for the basis function class. This is initialised with a scaler and
  parameter manager, and uses these to calculate the scale factors and
  derivatives for each reflection based on the model components."""

  # First initialise a scaling model, scaler and parameter manager.
  (test_reflections, test_experiments, params) = (
    generated_refl(), generated_single_exp(), generated_KB_param)
  assert len(test_experiments) == 1
  assert len(test_reflections) == 1
  experiments = create_scaling_model(params, test_experiments, test_reflections)
  scaler = create_scaler(params, experiments, test_reflections)
  assert scaler.experiments.scaling_model.id_ == 'KB'
  apm = scaling_active_parameter_manager(scaler.components, ['decay', 'scale'])

  # First test that scale factors can be successfully updated.
  # Manually change the parameters in the apm.
  n_scale_params = scaler.components['scale'].n_params
  n_decay_params = scaler.components['decay'].n_params
  # Note, order of params in apm.x depends on order in scaling model components.
  new_B = 1.0
  new_S = 2.0
  apm.x = (flex.double(n_scale_params, new_S))
  apm.x.extend(flex.double(n_decay_params, new_B))
  basis_fn = basis_function(scaler, apm)
  basis_fn.update_scale_factors_with_curvs()
  assert list(scaler.components['decay'].parameters) == (
    list(flex.double([new_B] * n_decay_params)))
  assert list(scaler.components['scale'].parameters) == (
    list(flex.double([new_S] * n_scale_params)))

  # Now test that the inverse scale factor is correctly calculated.
  calculated_sfs = basis_fn.calculate_scale_factors()
  assert list(calculated_sfs) == list(new_S * np.exp(new_B/
    (2.0*(scaler.components['decay'].d_values**2))))

  # Now check that the derivative matrix is correctly calculated.
  calc_derivs = basis_fn.calculate_derivatives()
  decay = scaler.components['decay']
  scale = scaler.components['scale']
  assert calc_derivs[0, 0] == scale.derivatives[0, 0] * decay.inverse_scales[0]
  assert calc_derivs[1, 0] == scale.derivatives[1, 0] * decay.inverse_scales[1]
  assert calc_derivs[2, 0] == scale.derivatives[2, 0] * decay.inverse_scales[2]
  assert calc_derivs[0, 1] == decay.derivatives[0, 0] * scale.inverse_scales[0]
  assert calc_derivs[1, 1] == decay.derivatives[1, 0] * scale.inverse_scales[1]
  assert calc_derivs[2, 1] == decay.derivatives[2, 0] * scale.inverse_scales[2]

  # Test that the curvatures matrix is correctly composed.
  calc_curvs = basis_fn.calculate_curvatures()
  assert calc_curvs[0, 0] == scale.curvatures[0, 0] * decay.inverse_scales[0]
  assert calc_curvs[1, 0] == scale.curvatures[1, 0] * decay.inverse_scales[1]
  assert calc_curvs[2, 0] == scale.curvatures[2, 0] * decay.inverse_scales[2]
  assert calc_curvs[0, 1] == decay.curvatures[0, 0] * scale.inverse_scales[0]
  assert calc_curvs[1, 1] == decay.curvatures[1, 0] * scale.inverse_scales[1]
  assert calc_curvs[2, 1] == decay.curvatures[2, 0] * scale.inverse_scales[2]

  # Repeat the test when there is only one active parameter.
  # First reset the parameters
  scaler.components['decay'].parameters = flex.double([0.0])
  scaler.components['scale'].parameters = flex.double([1.0])
  scaler.components['decay'].calculate_scales_and_derivatives()
  scaler.components['scale'].calculate_scales_and_derivatives()

  # Now generate a parameter manager for a single component.
  apm = scaling_active_parameter_manager(scaler.components, ['scale'])
  new_S = 2.0
  apm.x = flex.double(scaler.components['scale'].n_params, new_S)
  basis_func = basis_function(scaler, apm)
  basis_fn = basis_func.return_basis() # All in one alternative call.

  # Test that the scales and derivatives were correctly calculated
  assert list(basis_fn[0]) == list([new_S] *
    scaler.components['scale'].inverse_scales.size())
  assert basis_fn[1][0, 0] == scaler.components['scale'].derivatives[0, 0]
  assert basis_fn[1][1, 0] == scaler.components['scale'].derivatives[1, 0]
  assert basis_fn[1][2, 0] == scaler.components['scale'].derivatives[2, 0]
