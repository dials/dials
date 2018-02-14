'''
This code tests the data managers and active parameter managers.
'''
import copy as copy
import numpy as np
import pytest
from dials.array_family import flex
from dials.util.options import OptionParser
from ParameterHandler import scaling_active_parameter_manager
from active_parameter_managers import (multi_active_parameter_manager,
  active_parameter_manager)
from basis_functions import basis_function
from libtbx import phil
from dxtbx.model.experiment_list import ExperimentList
from dxtbx.model import Crystal, Scan, Beam, Goniometer, Detector, Experiment
from dials.algorithms.scaling.model import ScalingModelFactory as ScalingModelFactory
import ScalerFactory as ScalerFactory


def generated_refl():
  '''function to generate input for datamanagers'''
  #these miller_idx/d_values don't make physical sense, but I didn't want to
  #have to write the tests for lots of reflections.
  reflections = flex.reflection_table()
  reflections['intensity.prf.value'] = flex.double([1.0, 10.0, 100.0])
  reflections['intensity.prf.variance'] = flex.double([1.0, 10.0, 100.0])
  reflections['miller_index'] = flex.miller_index([(1, 0, 0), (0, 0, 1), (1, 0, 0)]) #don't change
  reflections['d'] = flex.double([0.8, 2.0, 2.0]) #don't change
  reflections['lp'] = flex.double([1.0, 1.0, 1.0])
  reflections['dqe'] = flex.double([1.0, 1.0, 1.0])
  reflections['partiality'] = flex.double([1.0, 1.0, 1.0])
  reflections['xyzobs.px.value'] = flex.vec3_double([(0.0, 0.0, 0.0),
    (0.0, 0.0, 5.0), (0.0, 0.0, 10.0)])
  reflections['s1'] = flex.vec3_double([(0.0, 0.1, 1.0), (0.0, 0.1, 1.0),
    (0.0, 0.1, 1.0)])
  reflections.set_flags(flex.bool([True, True, True]), reflections.flags.integrated)
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
  detector = Detector()
  experiments.append(Experiment(beam=beam, scan=scan, goniometer=goniometer,
    detector=detector, crystal=crystal))
  experiments.append(Experiment(beam=beam, scan=scan, goniometer=goniometer,
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
  parameters.__inject__('scaling_model', 'aimless')
  return parameters

@pytest.fixture(scope='module')
def generated_KB_param():
  phil_scope = phil.parse('''
      include scope dials.algorithms.scaling.scaling_options.phil_scope
  ''', process_includes=True)

  optionparser = OptionParser(phil=phil_scope, check_format=False)
  parameters, _ = optionparser.parse_args(args=None, quick_parse=True,
    show_diff_phil=False)
  parameters.__inject__('scaling_model', 'KB')
  return parameters

#@pytest.fixture(scope='module')
def generated_target_input(generated_refl, generated_two_exp, generated_param):
  refl = generated_refl
  exp = generated_two_exp
  param = generated_param
  refl_2 = copy.deepcopy(generated_refl[0])
  refl_2['inverse_scale_factor'] = flex.double([1.0, 1.0, 1.0])
  refl.append(refl_2)
  return (refl, exp, param)

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
  experiments = ScalingModelFactory.Factory.create(params, test_experiments, test_reflections)
  scaler = ScalerFactory.Factory.create(params, experiments, test_reflections)

  assert 'asu_miller_index' in scaler.reflection_table.keys()
  assert 'intensity' in scaler.reflection_table.keys()

  #test for successful escape from normalised intensity calculation
  n = (scaler.reflection_table['Esq'] == 0.0).count(True)
  assert n != len(scaler.reflection_table)

def test_TargetScalerFactory():
  '''test for successful initialisation of targetedDataManager'''
  (test_reflections, test_experiments, params) = generated_target_input(
    generated_refl(), generated_single_exp(), generated_param())

  experiments = ScalingModelFactory.Factory.create(params, test_experiments, test_reflections)
  _ = ScalerFactory.Factory.create(params, experiments, test_reflections)

def test_apm():
  '''test for a single active parameter manager. Also serves as general
  test for initialisation of AimlessDataManager'''
  (test_reflections, test_experiments, params) = generated_single_input(
    generated_refl(), generated_single_exp(), generated_param())
  assert len(test_experiments) == 1
  assert len(test_reflections) == 1
  experiments = ScalingModelFactory.Factory.create(params, test_experiments, test_reflections)
  scaler = ScalerFactory.Factory.create(params, experiments, test_reflections)

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

  experiments = ScalingModelFactory.Factory.create(params, test_experiments, test_reflections)
  scaler = ScalerFactory.Factory.create(params, experiments, test_reflections)

  assert isinstance(scaler, ScalerFactory.MultiScaler)
  components_list = [i.components for i in scaler.single_scalers]
  apm = multi_active_parameter_manager(components_list,
    [['decay', 'scale'], ['decay', 'scale']], scaling_active_parameter_manager)
  assert 'decay' in apm.components_list
  assert 'scale' in apm.components_list
  assert 'absorption' not in apm.components_list
  dm0 = scaler.single_scalers[0]
  dm1 = scaler.single_scalers[1]
  assert apm.n_active_params == (dm0.components['scale'].n_params + dm0.components['decay'].n_params
    + dm1.components['scale'].n_params + dm1.components['decay'].n_params)
  assert apm.apm_data[0]['start_idx'] == 0
  assert apm.apm_data[0]['end_idx'] == apm.apm_list[0].n_active_params
  assert apm.apm_data[1]['start_idx'] == apm.apm_list[0].n_active_params
  assert apm.apm_data[1]['end_idx'] == apm.n_active_params
  params = apm.apm_list[0].x
  params.extend(apm.apm_list[1].x)
  assert list(apm.x) == list(params)


def test_basis_function(generated_KB_param):
  '''test functionality of basis function object'''
  (test_reflections, test_experiments, params) = generated_single_input(
    generated_refl(), generated_single_exp(), generated_KB_param)
  assert len(test_experiments) == 1
  assert len(test_reflections) == 1
  experiments = ScalingModelFactory.Factory.create(params, test_experiments, test_reflections)
  scaler = ScalerFactory.Factory.create(params, experiments, test_reflections)
  assert scaler.id_ == 'KB'

  #(test_reflections, test_experiments, params) = generated_input
  #params.scaling_options.multi_mode = False
  #data_manager = KB_Data_Manager(test_reflections, test_experiments, params)

  apm = scaling_active_parameter_manager(scaler.components, ['decay', 'scale'])

  #first change the parameters in the apm
  n_scale_params = scaler.components['scale'].n_params
  n_decay_params = scaler.components['decay'].n_params
  #order of params in apm.x - depends on order in scaling model 'components'
  new_B = 1.0
  new_S = 2.0
  apm.x = (flex.double([new_S] * n_scale_params))
  apm.x.extend(flex.double([new_B] * n_decay_params))
  basis_func = basis_function(scaler, apm)

  #first check that scale factors can be successfully updated
  basis_fn = basis_func.return_basis() #includes update SF method.
  assert list(scaler.components['decay'].parameters) == (
    list(flex.double([new_B] * n_decay_params)))
  assert list(scaler.components['scale'].parameters) == (
    list(flex.double([new_S] * n_scale_params)))

  #now check that the inverse scale factor was correctly calculated
  assert list(basis_fn[0]) == list(new_S * np.exp(new_B/
    (2.0*(scaler.components['decay'].d_values**2))))

  #now check that the derivative matrix was correctly calculated
  g_decay = scaler.components['decay']
  g_scale = scaler.components['scale']
  assert basis_fn[1][0, 0] == g_scale.derivatives[0, 0] * g_decay.inverse_scales[0]
  assert basis_fn[1][1, 0] == g_scale.derivatives[1, 0] * g_decay.inverse_scales[1]
  assert basis_fn[1][2, 0] == g_scale.derivatives[2, 0] * g_decay.inverse_scales[2]
  assert basis_fn[1][0, 1] == g_decay.derivatives[0, 0] * g_scale.inverse_scales[0]
  assert basis_fn[1][1, 1] == g_decay.derivatives[1, 0] * g_scale.inverse_scales[1]
  assert basis_fn[1][2, 1] == g_decay.derivatives[2, 0] * g_scale.inverse_scales[2]

  #test case of only one active parameter as well on fresh data manager
  (test_reflections, test_experiments, params) = generated_single_input(
    generated_refl(), generated_single_exp(), generated_KB_param)
  assert len(test_experiments) == 1
  assert len(test_reflections) == 1
  experiments = ScalingModelFactory.Factory.create(params, test_experiments, test_reflections)
  scaler = ScalerFactory.Factory.create(params, experiments, test_reflections)
  apm = scaling_active_parameter_manager(scaler.components, ['scale'])

  #need to set the inverse scales for this test dataset.
  scaler.components['scale'].inverse_scales = flex.double([1.0, 1.0, 1.0])
  scaler.components['decay'].inverse_scales = flex.double([1.0, 1.0, 1.0])
  #first change the parameters in the apm
  new_S = 2.0
  apm.x = flex.double([new_S] * scaler.components['scale'].n_params)
  basis_func = basis_function(scaler, apm)
  basis_fn = basis_func.return_basis()

  #now check that the inverse scale factor was correctly calculated
  assert list(basis_fn[0]) == list([new_S] * len(scaler.components['scale'].inverse_scales))

  #now check that the derivative matrix was correctly calculated
  assert basis_fn[1][0, 0] == scaler.components['scale'].derivatives[0, 0]
  assert basis_fn[1][1, 0] == scaler.components['scale'].derivatives[1, 0]
  assert basis_fn[1][2, 0] == scaler.components['scale'].derivatives[2, 0]

def test_target_function(generated_KB_param):
  (test_reflections, test_experiments, params) = generated_single_input(
    generated_refl(), generated_single_exp(), generated_KB_param)
  assert len(test_experiments) == 1
  assert len(test_reflections) == 1
  experiments = ScalingModelFactory.Factory.create(params, test_experiments, test_reflections)
  scaler = ScalerFactory.Factory.create(params, experiments, test_reflections)
  assert scaler.id_ == 'KB'

  #setup of data manager
  scaler.components['scale'].parameters = flex.double([2.0])
  scaler.components['decay'].parameters = flex.double([1.0])
  scaler.components['scale'].inverse_scales = flex.double([1.0, 2.0, 1.0])
  scaler.components['decay'].inverse_scales = flex.double([1.0, 2.0, 1.0])

  apm = scaling_active_parameter_manager(scaler.components, ['scale', 'decay'])
  scaler.update_for_minimisation(apm)

  #first get residual and gradients
  res, grad = scaler.get_target_function(apm)
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
