'''
This code tests the data managers and active parameter managers.
'''
import copy as copy
import numpy as np
import pytest
from dials.array_family import flex
from dials.util.options import OptionParser
from data_manager_functions import (ScalingDataManager, MultiCrystalDataManager,
  AimlessDataManager, TargetedDataManager, KB_Data_Manager,
  active_parameter_manager, multi_active_parameter_manager)
from basis_functions import basis_function
from target_function import target_function
from libtbx import phil
from dxtbx.model.experiment_list import ExperimentList
from dxtbx.model import Crystal, Scan, Beam, Goniometer

@pytest.fixture(scope='module')
def generated_input():
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

  experiments = ExperimentList()
  exp_dict = {"__id__" : "crystal", "real_space_a": [1.0, 0.0, 0.0],
              "real_space_b": [0.0, 1.0, 0.0], "real_space_c": [0.0, 0.0, 2.0],
              "space_group_hall_symbol": " C 2y"}
  experiments.crystal = Crystal.from_dict(exp_dict)
  experiments.scan = Scan(image_range=[0, 90], oscillation=[0.0, 1.0])
  experiments.beam = Beam(s0=(0.0, 0.0, 1.01))
  experiments.goniometer = Goniometer((1.0, 0.0, 0.0))

  phil_scope = phil.parse('''
      include scope dials.algorithms.scaling.scaling_options.phil_scope
  ''', process_includes=True)

  optionparser = OptionParser(phil=phil_scope, check_format=False)
  parameters, _ = optionparser.parse_args(args=None, quick_parse=True,
    show_diff_phil=False)
  parameters.__inject__('scaling_method', 'aimless')
  return (reflections, experiments, parameters)


def test_ScalingDataManager(generated_input):
  '''test a few extra features that are not covered by test_Ih_table'''
  (test_reflections, test_experiments, params) = generated_input
  data_manager = ScalingDataManager(test_reflections, test_experiments, params)

  assert 'asu_miller_index' in data_manager.reflection_table.keys()
  assert 'intensity' in data_manager.reflection_table.keys()

  #test for successful escape from normalised intensity calculation
  n = (data_manager.reflection_table['Esq'] == 0.0).count(True)
  assert n != len(data_manager.reflection_table)

'''def test_KBDataManager():
  #test for successful initialisation of KBDataManager
  (test_reflections, test_experiments, params) = generate_test_input()
  data_manager = KB_Data_Manager(test_reflections, test_experiments, params)'''

'''def test_AimlessDataManager():
  (test_reflections, test_experiments, params) = generate_test_input()
  data_manager = AimlessDataManager(test_reflections, test_experiments, params)'''


def test_targeted_data_manager(generated_input):
  '''test for successful initialisation of targetedDataManager'''
  (test_reflections, test_experiments, params) = generated_input
  targeted_reflections = copy.deepcopy(test_reflections)
  targeted_reflections['inverse_scale_factor'] = flex.double([1.0, 1.0, 1.0])
  targeted_dm = TargetedDataManager(test_reflections, test_experiments,
    targeted_reflections, params)

  # note - write tests to check correct behaviour w.r.t Ih_table setup here.
  #also do a test update cycle and check Ih values haven't been overwritten.

'''def test_multi_data_manager():
  (test_reflections, test_experiments, params) = generate_test_input()
  params.scaling_options.__inject__('multi_mode', True)
  data_manager = MultiCrystalDataManager([test_reflections, test_reflections],
    [test_experiments, test_experiments], params)'''

def test_apm(generated_input):
  '''test for a single active parameter manager. Also serves as general
  test for initialisation of AimlessDataManager'''
  (test_reflections, test_experiments, params) = generated_input
  params.scaling_options.__inject__('multi_mode', False)
  data_manager = AimlessDataManager(test_reflections, test_experiments, params)

  apm = active_parameter_manager(data_manager, ['g_decay', 'g_scale'])
  assert 'g_decay' in apm.active_parameterisation
  assert 'g_scale' in apm.active_parameterisation
  assert 'g_absorption' not in apm.active_parameterisation
  assert apm.n_active_params == data_manager.g_scale.n_params + data_manager.g_decay.n_params
  assert apm.constant_g_values is not None
  assert list(apm.constant_g_values) == list(data_manager.g_absorption.inverse_scales)
  for i, p in enumerate(apm.active_parameterisation):
    assert apm.n_params_list[i] == data_manager.g_parameterisation[p].n_params

def test_multi_apm(generated_input):
  '''test for the multi active parameter manager. Also serves as general
  test for initialisation of MultiCrystalDataManager'''
  (test_reflections, test_experiments, params) = generated_input
  params.scaling_options.multi_mode = True
  data_manager = MultiCrystalDataManager([test_reflections, test_reflections],
    [test_experiments, test_experiments], params)

  apm = multi_active_parameter_manager(data_manager, ['g_decay', 'g_scale'])
  assert 'g_decay' in apm.active_parameterisation
  assert 'g_scale' in apm.active_parameterisation
  assert 'g_absorption' not in apm.active_parameterisation
  dm0 = data_manager.data_managers[0]
  dm1 = data_manager.data_managers[1]
  assert apm.n_active_params == (dm0.g_scale.n_params + dm0.g_decay.n_params
    + dm1.g_scale.n_params + dm1.g_decay.n_params)
  assert apm.n_params_in_each_apm[0] == dm0.g_scale.n_params + dm0.g_decay.n_params
  assert apm.n_params_in_each_apm[1] == dm1.g_scale.n_params + dm1.g_decay.n_params
  params = apm.apm_list[0].x
  params.extend(apm.apm_list[1].x)
  assert list(apm.x) == list(params)

def test_basis_function(generated_input):
  '''test functionality of basis function object'''
  (test_reflections, test_experiments, params) = generated_input
  params.scaling_options.multi_mode = False
  data_manager = KB_Data_Manager(test_reflections, test_experiments, params)

  apm = active_parameter_manager(data_manager, ['g_scale', 'g_decay'])

  #first change the parameters in the apm
  n_scale_params = data_manager.g_scale.n_params
  n_decay_params = data_manager.g_decay.n_params
  #order of params in apm.x - decay first then scale
  new_B = 1.0
  new_S = 2.0
  apm.x = flex.double([new_B] * n_decay_params)
  apm.x.extend(flex.double([new_S] * n_scale_params))
  basis_func = basis_function(data_manager, apm)

  #first check that scale factors can be successfully updated
  basis_fn = basis_func.return_basis() #includes update SF method.
  assert list(data_manager.g_decay.parameters) == (
    list(flex.double([new_B] * n_decay_params)))
  assert list(data_manager.g_scale.parameters) == (
    list(flex.double([new_S] * n_scale_params)))

  #now check that the inverse scale factor was correctly calculated
  assert list(basis_fn[0]) == list(new_S * np.exp(new_B/
    (2.0*(data_manager.g_decay.d_values**2))))

  #now check that the derivative matrix was correctly calculated
  g_decay = data_manager.g_decay
  g_scale = data_manager.g_scale
  assert basis_fn[1][0, 0] == g_decay.derivatives[0, 0] * g_scale.inverse_scales[0]
  assert basis_fn[1][1, 0] == g_decay.derivatives[1, 0] * g_scale.inverse_scales[1]
  assert basis_fn[1][2, 0] == g_decay.derivatives[2, 0] * g_scale.inverse_scales[2]
  assert basis_fn[1][0, 1] == g_scale.derivatives[0, 0] * g_decay.inverse_scales[0]
  assert basis_fn[1][1, 1] == g_scale.derivatives[1, 0] * g_decay.inverse_scales[1]
  assert basis_fn[1][2, 1] == g_scale.derivatives[2, 0] * g_decay.inverse_scales[2]

  #test case of only one active parameter as well on fresh data manager
  data_manager = KB_Data_Manager(test_reflections, test_experiments, params)
  apm = active_parameter_manager(data_manager, ['g_scale'])

  #need to set the inverse scales for this test dataset.
  data_manager.g_scale.inverse_scales = flex.double([1.0, 1.0, 1.0])
  data_manager.g_decay.inverse_scales = flex.double([1.0, 1.0, 1.0])
  #first change the parameters in the apm
  new_S = 2.0
  apm.x = flex.double([new_S] * data_manager.g_scale.n_params)
  basis_func = basis_function(data_manager, apm)
  basis_fn = basis_func.return_basis()

  #now check that the inverse scale factor was correctly calculated
  assert list(basis_fn[0]) == list([new_S] * len(data_manager.g_scale.inverse_scales))

  #now check that the derivative matrix was correctly calculated
  assert basis_fn[1][0, 0] == data_manager.g_scale.derivatives[0, 0]
  assert basis_fn[1][1, 0] == data_manager.g_scale.derivatives[1, 0]
  assert basis_fn[1][2, 0] == data_manager.g_scale.derivatives[2, 0]

def test_target_function(generated_input):
  (test_reflections, test_experiments, params) = generated_input
  #setup of data manager
  data_manager = KB_Data_Manager(test_reflections, test_experiments, params)
  data_manager.g_scale.parameters = flex.double([2.0])
  data_manager.g_decay.parameters = flex.double([1.0])
  data_manager.g_scale.inverse_scales = flex.double([1.0, 2.0, 1.0])
  data_manager.g_decay.inverse_scales = flex.double([1.0, 2.0, 1.0])

  apm = active_parameter_manager(data_manager, ['g_scale', 'g_decay'])
  data_manager.update_for_minimisation(apm)

  #first get residual and gradients
  res, grad = data_manager.get_target_function(apm)
  assert res > 1e-8, """residual should not be zero, or the gradient test
    below will not really be working!"""
  f_d_grad = calculate_gradient_fd(data_manager, apm)
  assert list(grad - f_d_grad) < ([1e-5] * apm.n_active_params)
  assert list(f_d_grad) > [1e-8] * apm.n_active_params, """finite difference
    grad should not be zero, or the test will not really be working!"""

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
