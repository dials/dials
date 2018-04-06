'''
This code tests the data managers and active parameter managers.
'''
import copy as copy
import pytest
from scitbx import sparse
from dials.array_family import flex
from dials.util.options import OptionParser
from parameter_handler import scaling_active_parameter_manager
from active_parameter_managers import (multi_active_parameter_manager,
  active_parameter_manager)
from libtbx import phil
from libtbx.test_utils import approx_equal
from dxtbx.model.experiment_list import ExperimentList
from dxtbx.model import Crystal, Scan, Beam, Goniometer, Detector, Experiment
from dials.algorithms.scaling.model.scaling_model_factory import \
  create_scaling_model
from dials.algorithms.scaling.scaler_factory import create_scaler,\
  MultiScaler
from dials.algorithms.scaling.scaler import SingleScalerBase
from dials.algorithms.scaling.target_function import ScalingTarget,\
  ScalingTargetFixedIH

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




def test_apm():
  """test for a single active parameter manager. Also serves as general
  test for initialisation of PhysicalScaler."""
  (test_reflections, test_experiments, params) = generated_single_input(
    generated_refl(), generated_single_exp(), generated_param())
  assert len(test_experiments) == 1
  assert len(test_reflections) == 1
  experiments = create_scaling_model(params, test_experiments, test_reflections)
  scaler = create_scaler(params, experiments, test_reflections)

  apm = scaling_active_parameter_manager(scaler.components, ['decay', 'scale'])
  assert 'decay' in apm.components_list
  assert 'scale' in apm.components_list
  assert 'absorption' not in apm.components_list
  assert apm.n_active_params == (scaler.components['scale'].n_params
    + scaler.components['decay'].n_params)
  assert apm.constant_g_values is not None
  assert list(apm.constant_g_values) == list(scaler.components['absorption'].inverse_scales)
  for component in apm.components:
    assert apm.components[component]['n_params'] == scaler.components[component].n_params

def test_multi_apm():
  '''test for the multi active parameter manager. Also serves as general
  test for initialisation of MultiCrystalDataManager'''
  (test_reflections, test_experiments, params) = generated_multi_input(
    generated_refl(), generated_two_exp(), generated_param())

  experiments = create_scaling_model(params, test_experiments, test_reflections)
  scaler = create_scaler(params, experiments, test_reflections)

  assert isinstance(scaler, MultiScaler)
  components_list = [i.components for i in scaler.single_scalers]
  apm = multi_active_parameter_manager(components_list,
    [['decay', 'scale'], ['decay', 'scale']], scaling_active_parameter_manager)
  assert 'decay' in apm.components_list
  assert 'scale' in apm.components_list
  assert 'absorption' not in apm.components_list
  dm0 = scaler.single_scalers[0]
  dm1 = scaler.single_scalers[1]
  assert apm.n_active_params == (dm0.components['scale'].n_params
    + dm0.components['decay'].n_params + dm1.components['scale'].n_params
    + dm1.components['decay'].n_params)
  assert apm.apm_data[0]['start_idx'] == 0
  assert apm.apm_data[0]['end_idx'] == apm.apm_list[0].n_active_params
  assert apm.apm_data[1]['start_idx'] == apm.apm_list[0].n_active_params
  assert apm.apm_data[1]['end_idx'] == apm.n_active_params
  params = apm.apm_list[0].x
  params.extend(apm.apm_list[1].x)
  assert list(apm.x) == list(params)

def test_target_function(generated_KB_param):
  (test_reflections, test_experiments, params) = generated_single_input(
    generated_refl(), generated_single_exp(), generated_KB_param)
  assert len(test_experiments) == 1
  assert len(test_reflections) == 1
  experiments = create_scaling_model(params, test_experiments, test_reflections)
  scaler = create_scaler(params, experiments, test_reflections)
  assert scaler.experiments.scaling_model.id_ == 'KB'

  #setup of data manager
  scaler.components['scale'].parameters = flex.double([2.0])
  scaler.components['decay'].parameters = flex.double([1.0])
  scaler.components['scale'].inverse_scales = flex.double([1.0, 2.0, 1.0])
  scaler.components['decay'].inverse_scales = flex.double([1.0, 2.0, 1.0])

  apm = scaling_active_parameter_manager(scaler.components, ['scale', 'decay'])
  scaler.update_for_minimisation(apm)

  #first get residual and gradients
  target = ScalingTarget(scaler, apm)
  res, grad = target.compute_functional_gradients()
  assert res > 1e-8, """residual should not be zero, or the gradient test
    below will not really be working!"""
  f_d_grad = calculate_gradient_fd(scaler, apm)
  assert list(grad - f_d_grad) < ([1e-5] * apm.n_active_params)
  sel = f_d_grad > 1e-8
  assert sel, """assert sel has some elements, as finite difference grad should
    not all be zero, or the test will not really be working!
    (expect one to be zero for KB scaling example?)"""

def calculate_gradient_fd(dm, apm):
  '''calculate gradient array with finite difference approach'''
  delta = 1.0e-6
  Ih_tab = dm.Ih_table
  gradients = flex.double([0.0] * apm.n_active_params)
  #iterate over parameters, varying one at a time and calculating the gradient
  for i in range(apm.n_active_params):
    apm.x[i] -= 0.5 * delta
    dm.update_for_minimisation(apm)
    R_low = ((((Ih_tab.intensities - (Ih_tab.inverse_scale_factors
      * Ih_tab.Ih_values))**2) * Ih_tab.weights))
    apm.x[i] += delta
    dm.update_for_minimisation(apm)
    R_upper = ((((Ih_tab.intensities - (Ih_tab.inverse_scale_factors
      * Ih_tab.Ih_values))**2) * Ih_tab.weights))
    apm.x[i] -= 0.5 * delta
    dm.update_for_minimisation(apm)
    gradients[i] = (flex.sum(R_upper) - flex.sum(R_low)) / delta
  return gradients

def test_general_apm():
  'test for general active parameter manager'

  class dummy_obj(object):
    '''class to create dummy object with parameters attribute'''
    def __init__(self, parameters):
      self.parameters = parameters

  a = dummy_obj(flex.double([1.0, 2.0]))
  b = dummy_obj(flex.double([2.0, 4.0]))
  c = dummy_obj(flex.double([3.0, 6.0]))

  components = {'1' : a, '2' : b, '3' : c}

  apm = active_parameter_manager(components, ['1', '3'])
  assert apm.n_active_params == 4
  assert '1' in apm.components_list
  assert '3' in apm.components_list
  assert '2' not in apm.components_list
  assert '1' in apm.components
  assert '3' in apm.components
  assert '2' not in apm.components
  assert apm.components['1']['object'].parameters == flex.double([1.0, 2.0])
  assert apm.components['1']['object'].parameters == flex.double([3.0, 6.0])

  assert apm.select_parameters('1') == flex.double([1.0, 2.0])
  assert apm.select_parameters('3') == flex.double([3.0, 6.0])

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
  assert list(components['scale'].derivatives.col(0)) == [(0,1.0), (1,1.0), (2,1.0)]
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
