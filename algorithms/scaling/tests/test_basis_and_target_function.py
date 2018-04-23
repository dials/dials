"""
Test for the basis function and target function module.
"""
import copy
import numpy as np
import pytest
from scitbx import sparse
from dials.array_family import flex
from dials.util.options import OptionParser
from parameter_handler import scaling_active_parameter_manager
from libtbx import phil
from libtbx.test_utils import approx_equal
from dxtbx.model.experiment_list import ExperimentList
from dxtbx.model import Crystal, Scan, Beam, Goniometer, Detector, Experiment
from dials.algorithms.scaling.scaling_library import create_scaling_model
from dials.algorithms.scaling.scaler_factory import create_scaler
from dials.algorithms.scaling.target_function import ScalingTarget
from dials.algorithms.scaling.basis_functions import basis_function
from dials.algorithms.scaling.model.components.scale_components import \
  SingleBScaleFactor, SingleScaleFactor

@pytest.fixture
def large_reflection_table():
  """Create a larger reflection table"""
  return generated_10_refl()

@pytest.fixture
def small_reflection_table():
  """Create a small reflection table"""
  return generated_refl()

def generated_10_refl():
  """Generate reflection table to test the basis and target function."""
  #these miller_idx/d_values don't make physical sense, but I didn't want to
  #have to write the tests for lots of reflections.
  reflections = flex.reflection_table()
  reflections['intensity.prf.value'] = flex.double([75.0, 10.0, 100.0, 25.0, 50.0, 100.0,
    25.0, 20.0, 300.0, 10.0])
  reflections['intensity.prf.variance'] = flex.double([50.0, 10.0, 100.0, 50.0, 10.0, 100.0,
    50.0, 10.0, 100.0, 10.0])
  reflections['miller_index'] = flex.miller_index([(1, 0, 0), (0, 0, 1),
    (1, 0, 0), (1, 0, 0), (0, 0, 1),
    (1, 0, 0), (0, 4, 0), (0, 0, 1),
    (1, 0, 0), (0, 4, 0)]) #don't change
  reflections['d'] = flex.double([2.0, 0.8, 2.0, 2.0, 0.8, 2.0, 2.0, 0.8, 2.0, 1.0]) #don't change
  reflections['lp'] = flex.double(10, 1.0)
  reflections['dqe'] = flex.double(10, 1.0)
  reflections['partiality'] = flex.double(10, 1.0)
  reflections['xyzobs.px.value'] = flex.vec3_double([(0.0, 0.0, 0.0),
    (0.0, 0.0, 5.0), (0.0, 0.0, 10.0), (0.0, 0.0, 15.0), (0.0, 0.0, 20.0),
    (0.0, 0.0, 25.0), (0.0, 0.0, 30.0), (0.0, 0.0, 35.0), (0.0, 0.0, 40.0),
    (0.0, 0.0, 59.0)])
  reflections['s1'] = flex.vec3_double([(0.0, 0.1, 1.0), (0.0, 0.1, 1.0),
    (0.0, 0.1, 20.0), (0.0, 0.1, 20.0), (0.0, 0.1, 20.0), (0.0, 0.1, 20.0),
    (0.0, 0.1, 20.0), (0.0, 0.1, 20.0), (0.0, 0.1, 20.0), (0.0, 0.1, 20.0)])
  reflections.set_flags(flex.bool(10, True), reflections.flags.integrated)
  return [reflections]


def generated_refl():
  """Generate reflection table to test the basis and target function."""
  #these miller_idx/d_values don't make physical sense, but I didn't want to
  #have to write the tests for lots of reflections.
  reflections = flex.reflection_table()
  reflections['intensity.prf.value'] = flex.double([75.0, 10.0, 100.0])
  reflections['intensity.prf.variance'] = flex.double([50.0, 10.0, 100.0])
  reflections['miller_index'] = flex.miller_index([(1, 0, 0), (0, 0, 1),
    (1, 0, 0)]) #don't change
  reflections['d'] = flex.double([2.0, 0.8, 2.0]) #don't change
  reflections['lp'] = flex.double([1.0, 1.0, 1.0])
  reflections['dqe'] = flex.double([1.0, 1.0, 1.0])
  reflections['partiality'] = flex.double([1.0, 1.0, 1.0])
  reflections['xyzobs.px.value'] = flex.vec3_double([(0.0, 0.0, 0.0),
    (0.0, 0.0, 5.0), (0.0, 0.0, 10.0)])
  reflections['s1'] = flex.vec3_double([(0.0, 0.1, 1.0), (0.0, 0.1, 1.0),
    (0.0, 0.1, 20.0)])
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
  scan = Scan(image_range=[0, 60], oscillation=[0.0, 1.0])
  beam = Beam(s0=(0.0, 0.0, 1.01))
  goniometer = Goniometer((1.0, 0.0, 0.0))
  detector = Detector()
  experiments.append(Experiment(beam=beam, scan=scan, goniometer=goniometer,
    detector=detector, crystal=crystal))
  return experiments

def generated_param(model='KB'):
  """Generate the scaling phil param scope."""
  phil_scope = phil.parse('''
      include scope dials.algorithms.scaling.scaling_options.phil_scope
  ''', process_includes=True)

  optionparser = OptionParser(phil=phil_scope, check_format=False)
  parameters, _ = optionparser.parse_args(args=None, quick_parse=True,
    show_diff_phil=False)
  parameters.__inject__('model', model)
  parameters.parameterisation.absorption_term = False
  return parameters

@pytest.fixture
def jacobian_gradient_input(large_reflection_table):
  return generated_param(model='physical'), generated_single_exp(), large_reflection_table


def test_basis_function(small_reflection_table):
  """Test for the basis function class. This calculates scale factors and
  derivatives for reflections based on the model components."""

  # To test the basis function, need a scaling active parameter manager - to set
  # this up we need a components dictionary with some reflection data.

  # Let's use KB model components for simplicity.
  rt = small_reflection_table[0]
  components = {'scale' : SingleScaleFactor(flex.double([1.0])), 'decay':
    SingleBScaleFactor(flex.double([0.0]))} #Create empty components.
  for component in components.itervalues():
    component.update_reflection_data(rt) #Add some data to components.

  apm = scaling_active_parameter_manager(components, ['decay', 'scale'])

  # First test that scale factors can be successfully updated.
  # Manually change the parameters in the apm.
  decay = components['decay'] # Define alias
  scale = components['scale'] # Define alias
  # Note, order of params in apm.x depends on order in scaling model components.
  new_B = 1.0
  new_S = 2.0
  apm.set_param_vals(flex.double([new_S, new_B]))
  basis_fn = basis_function(apm, curvatures=True)
  basis_fn.update_scale_factors()

  # Now test that the inverse scale factor is correctly calculated.
  calculated_sfs = basis_fn.calculate_scale_factors()
  assert list(calculated_sfs) == list(new_S * np.exp(new_B/
    (2.0*(decay.d_values**2))))

  # Now check that the derivative matrix is correctly calculated.
  calc_derivs = basis_fn.calculate_derivatives()
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
  components['decay'].parameters = flex.double([0.0])
  components['scale'].parameters = flex.double([1.0])
  components['decay'].calculate_scales_and_derivatives()
  components['scale'].calculate_scales_and_derivatives()

  # Now generate a parameter manager for a single component.
  apm = scaling_active_parameter_manager(components, ['scale'])
  new_S = 2.0
  apm.set_param_vals(flex.double(components['scale'].n_params, new_S))
  basis_fn = basis_function(apm).return_basis()
  #basis_fn = basis_func.return_basis() # All in one alternative call.

  # Test that the scales and derivatives were correctly calculated
  assert list(basis_fn[0]) == list([new_S] *
    components['scale'].inverse_scales.size())
  assert basis_fn[1][0, 0] == components['scale'].derivatives[0, 0]
  assert basis_fn[1][1, 0] == components['scale'].derivatives[1, 0]
  assert basis_fn[1][2, 0] == components['scale'].derivatives[2, 0]

def test_target_function():
  """Test for the ScalingTarget class."""

  # First set up the scaler
  (test_reflections, test_experiments, params) = (
    generated_refl(), generated_single_exp(), generated_param(model='KB'))
  assert len(test_experiments) == 1
  assert len(test_reflections) == 1
  experiments = create_scaling_model(params, test_experiments, test_reflections)
  scaler = create_scaler(params, experiments, test_reflections)
  assert scaler.experiments.scaling_model.id_ == 'KB'

  # Initialise the parameters and create an apm
  scaler.components['scale'].parameters = flex.double([2.0])
  scaler.components['decay'].parameters = flex.double([0.0])
  scaler.components['scale'].inverse_scales = flex.double([2.0, 2.0, 2.0])
  scaler.components['decay'].inverse_scales = flex.double([1.0, 1.0, 1.0])
  apm = scaling_active_parameter_manager(scaler.components, ['scale', 'decay'])
  scaler.update_for_minimisation(apm)

  # Create a scaling target and check gradients
  target = ScalingTarget(scaler, apm)
  res, grad = target.compute_functional_gradients()
  assert res > 1e-8, """residual should not be zero, or the gradient test
    below will not really be working!"""
  f_d_grad = calculate_gradient_fd(target)
  assert approx_equal(list(grad), list(f_d_grad))
  sel = f_d_grad > 1e-8
  assert sel, """assert sel has some elements, as finite difference grad should
    not all be zero, or the test will not really be working!
    (expect one to be zero for KB scaling example?)"""

  assert target.get_num_matches() == 3
  # Below methods needed for refinement engine calls
  _ = target.compute_restraints_residuals_and_gradients()
  _ = target.compute_residuals_and_gradients()
  _ = target.compute_residuals()
  _ = target.compute_functional_gradients()
  _ = target.achieved()
  _ = target.predict()
  resid = (target.calculate_residuals()**2) * target.weights
  # Note - activate two below when curvatures are implemented.
  #_ = target.compute_restraints_functional_gradients_and_curvatures()
  #_ = target.compute_functional_gradients_and_curvatures()

  # Calculate residuals explicitly and check RMSDS.
  assert approx_equal(list(resid), [50.0/36.0, 0.0, 100.0/36.0])
  assert approx_equal(target.rmsds()[0], (150.0/(36.0*3.0))**0.5)

def test_target_jacobian_calc(jacobian_gradient_input):
  test_params, exp, test_refl = jacobian_gradient_input
  test_params.parameterisation.decay_term=False
  experiments = create_scaling_model(test_params, exp, test_refl)
  assert experiments[0].scaling_model.id_ == 'physical'
  scaler = create_scaler(test_params, experiments, test_refl)

  apm = scaling_active_parameter_manager(scaler.components, ['scale'])

  target = ScalingTarget(scaler, apm)
  target.predict()

  fd_jacobian = calculate_jacobian_fd(target)
  print(fd_jacobian)
  r, jacobian, w = target.compute_residuals_and_gradients()
  r = (r/w)**0.5
  print(jacobian)
  print(list(w))
  for i in range(0, 3):
    for j in range(0, 2):
      assert approx_equal(jacobian[i, j], fd_jacobian[i, j])


def calculate_gradient_fd(target):
  """Calculate gradient array with finite difference approach."""
  delta = 1.0e-6
  gradients = flex.double([0.0] * target.apm.n_active_params)
  #iterate over parameters, varying one at a time and calculating the gradient
  for i in range(target.apm.n_active_params):
    new_x = copy.copy(target.apm.x)
    new_x[i] -= 0.5 * delta
    #target.apm.x[i] -= 0.5 * delta
    target.apm.set_param_vals(new_x)
    target.predict()
    R_low = (target.calculate_residuals()**2) * target.weights
    #target.apm.x[i] += delta
    new_x[i] += delta
    target.apm.set_param_vals(new_x)
    target.predict()
    R_upper = (target.calculate_residuals()**2) * target.weights
    #target.apm.x[i] -= 0.5 * delta
    new_x[i] -= 0.5 * delta
    target.apm.set_param_vals(new_x)
    target.predict()
    gradients[i] = (flex.sum(R_upper) - flex.sum(R_low)) / delta
  return gradients

def calculate_jacobian_fd(target):
  """Calculate jacobian matrix with finite difference approach."""
  delta = 1.0e-8
  #apm = target.apm
  jacobian = sparse.matrix(target.get_num_matches(), target.apm.n_active_params)
  #iterate over parameters, varying one at a time and calculating the residuals
  for i in range(target.apm.n_active_params):
    new_x = copy.copy(target.apm.x)
    new_x[i] -= 0.5 * delta
    target.apm.set_param_vals(new_x)
    target.predict()
    R_low = target.calculate_residuals()#unweighted unsquared residual
    new_x[i] += delta
    target.apm.set_param_vals(new_x)
    target.predict()
    R_upper = target.calculate_residuals() #unweighted unsquared residual
    new_x[i] -= 0.5 * delta
    target.apm.set_param_vals(new_x)
    target.predict()
    fin_difference = (R_upper - R_low) / delta
    for j in range(fin_difference.size()):
      jacobian[j, i] = fin_difference[j]
  return jacobian

