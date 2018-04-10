'''
This code tests the data managers and active parameter managers.
'''
import copy as copy
import pytest
from scitbx import sparse
from dials.array_family import flex
from dials.util.options import OptionParser
from libtbx import phil
from libtbx.test_utils import approx_equal
from dxtbx.model.experiment_list import ExperimentList
from dxtbx.model import Crystal, Scan, Beam, Goniometer, Detector, Experiment
from dials.algorithms.scaling.model.scaling_model_factory import \
  create_scaling_model
from dials.algorithms.scaling.scaler_factory import create_scaler
from dials.algorithms.scaling.target_function import ScalingTarget
from dials_scaling.jbe.scaling_code.scaler import SingleScalerBase

def generated_refl():
  '''function to generate input for datamanagers'''
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

#@pytest.fixture(scope='module')
def generated_two_exp():
  experiments = ExperimentList()
  exp_dict = {"__id__" : "crystal", "real_space_a": [1.0, 0.0, 0.0],
              "real_space_b": [0.0, 1.0, 0.0], "real_space_c": [0.0, 0.0, 2.0],
              "space_group_hall_symbol": " C 2y"}
  crystal = Crystal.from_dict(exp_dict)
  scan = Scan(image_range=[0, 90], oscillation=[0.0, 1.0])
  beam = Beam(s0=(0.0, 0.0, 1.01))
  goniometer = Goniometer((1.0, 0.0, 0.0))
  goniometer_2 = Goniometer((1.0, 1.0, 0.0))
  detector = Detector()
  experiments.append(Experiment(beam=beam, scan=scan, goniometer=goniometer,
    detector=detector, crystal=crystal))
  experiments.append(Experiment(beam=beam, scan=scan, goniometer=goniometer_2,
    detector=detector, crystal=crystal))
  return experiments

@pytest.fixture(scope='module')
def generated_param():
  phil_scope = phil.parse('''
      include scope dials.algorithms.scaling.scaling_options.phil_scope
  ''', process_includes=True)

  optionparser = OptionParser(phil=phil_scope, check_format=False)
  parameters, _ = optionparser.parse_args(args=None, quick_parse=True,
    show_diff_phil=False)
  parameters.__inject__('model', 'physical')
  return parameters

@pytest.fixture(scope='module')
def generated_KB_param():
  phil_scope = phil.parse('''
      include scope dials.algorithms.scaling.scaling_options.phil_scope
  ''', process_includes=True)

  optionparser = OptionParser(phil=phil_scope, check_format=False)
  parameters, _ = optionparser.parse_args(args=None, quick_parse=True,
    show_diff_phil=False)
  parameters.__inject__('model', 'KB')
  return parameters

#@pytest.fixture(scope='module')


def generated_multi_input(generated_refl, generated_two_exp, generated_param):
  refl = generated_refl
  exp = generated_two_exp
  param = generated_param
  refl_2 = copy.deepcopy(generated_refl[0])
  refl.append(refl_2)
  return (refl, exp, param)

#@pytest.fixture(scope='module')
def generated_single_input(generated_refl, generated_single_exp, generated_param):
  return (generated_refl, generated_single_exp, generated_param)

def test_SingleScalerFactory():
  '''test a few extra features that are not covered by test_Ih_table'''
  (test_reflections, test_experiments, params) = generated_single_input(
    generated_refl(), generated_single_exp(), generated_param())
  assert len(test_experiments) == 1
  assert len(test_reflections) == 1
  experiments = create_scaling_model(params, test_experiments, test_reflections)
  scaler = create_scaler(params, experiments, test_reflections)

  assert 'asu_miller_index' in scaler.reflection_table.keys()
  assert 'intensity' in scaler.reflection_table.keys()

  #test for successful escape from normalised intensity calculation
  n = (scaler.reflection_table['Esq'] == 0.0).count(True)
  assert n != len(scaler.reflection_table)

'''def calculate_jacobian_fd(target):
  """Calculate jacobian matrix with finite difference approach."""
  delta = 1.0e-6
  #apm = target.apm
  jacobian = sparse.matrix(target.get_num_matches(), target.apm.n_active_params)
  #iterate over parameters, varying one at a time and calculating the residuals
  for i in range(target.apm.n_active_params):
    target.apm.x[i] -= 0.5 * delta
    target.predict()
    R_low = (target.calculate_residuals()/target.weights)**0.5 #unweighted unsquared residual
    target.apm.x[i] += delta
    target.predict()
    R_upper = (target.calculate_residuals()/target.weights)**0.5 #unweighted unsquared residual
    target.apm.x[i] -= 0.5 * delta
    target.predict()
    fin_difference = (R_upper - R_low) / delta
    for j in range(fin_difference.size()):
      jacobian[j, i] = fin_difference[j]
  return jacobian

def test_target_jacobian_calc():
  (test_reflections, test_experiments, params) = generated_single_input(
    generated_refl(), generated_single_exp(), generated_param())
  assert len(test_experiments) == 1
  assert len(test_reflections) == 1
  experiments = create_scaling_model(params, test_experiments, test_reflections)
  scaler = create_scaler(params, experiments, test_reflections)

  apm = scaling_active_parameter_manager(scaler.components, ['decay', 'scale'])

  target = ScalingTarget(scaler, apm)

  fd_jacobian = calculate_jacobian_fd(target)
  print(fd_jacobian)
  r, jacobian, w = target.compute_residuals_and_gradients()
  print(jacobian)
  print(list(w))
  for i in range(0, 3):
    for j in range(0, 2):
      assert approx_equal(jacobian[i, j], fd_jacobian[i, j])'''

def test_sf_variance_calculation(generated_KB_param):
  (test_reflections, test_experiments, params) = generated_single_input(
    generated_refl(), generated_single_exp(), generated_KB_param)
  assert len(test_experiments) == 1
  assert len(test_reflections) == 1
  experiments = create_scaling_model(params, test_experiments, test_reflections)
  from dials.algorithms.scaling.scaler import calc_sf_variances
  components = experiments[0].scaling_model.components
  rt = flex.reflection_table()
  d1 = 1.0
  d2 = 2.0
  d3 = 3.0
  rt['d'] = flex.double([d1, d2, d3])
  components['scale'].update_reflection_data(rt)
  components['scale'].calculate_scales_and_derivatives()
  assert list(components['scale'].derivatives.col(0)) == [
    (0, 1.0), (1, 1.0), (2, 1.0)]
  components['decay'].update_reflection_data(rt)
  components['decay'].calculate_scales_and_derivatives()
  assert list(components['decay'].derivatives.col(0)) == [(0, 1.0/(2.0*d1*d1)),
    (1, 1.0/(2.0*d2*d2)), (2, 1.0/(2.0*d3*d3))]
  var_cov = sparse.matrix(2, 2)
  a = 0.2
  b = 0.3
  c = 0.1
  var_cov[0, 0] = a
  var_cov[0, 1] = c
  var_cov[1, 0] = c
  var_cov[1, 1] = b
  variances = calc_sf_variances(components, var_cov)
  assert approx_equal(list(variances), [b/(4.0*(d1**4.0)) + c/((d1**2.0)) + a,
    b/(4.0*(d2**4.0)) + c/((d2**2.0)) + a,
    b/(4.0*(d3**4.0)) + c/((d3**2.0)) + a])
