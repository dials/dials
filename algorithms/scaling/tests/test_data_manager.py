import copy as copy
from dials.array_family import flex
from dials.util.options import OptionParser
from data_manager_functions import (ScalingDataManager, MultiCrystalDataManager,
  AimlessDataManager, TargetedDataManager, KB_Data_Manager,
  active_parameter_manager)
from libtbx import phil
from dxtbx.model.experiment_list import ExperimentList
from dxtbx.model import Crystal, Scan, Beam, Goniometer

def generate_test_input():
  reflections = flex.reflection_table()
  reflections['intensity.prf.value'] = flex.double([1.0, 10.0, 100.0])
  reflections['intensity.prf.variance'] = flex.double([1.0, 10.0, 100.0])
  reflections['miller_index'] = flex.miller_index([(4, 0, 0), (0, 1, 0), (1, 0, 0)])
  reflections['d'] = flex.double([0.8, 2.0, 2.0])
  reflections['lp'] = flex.double([1.0, 1.0, 1.0])
  reflections['dqe'] = flex.double([1.0, 1.0, 1.0])
  reflections['partiality'] = flex.double([1.0, 1.0, 1.0])
  reflections['xyzobs.px.value'] = flex.vec3_double([(0.0, 0.0, 0.0),
    (0.0, 0.0, 45.0), (0.0, 0.0, 90.0)])
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


def test_ScalingDataManager():
  '''test a few extra features that are not covered by test_Ih_table'''
  (test_reflections, test_experiments, params) = generate_test_input()
  data_manager = ScalingDataManager(test_reflections, test_experiments, params)

  assert 'asu_miller_index' in data_manager.reflection_table.keys()
  assert 'intensity' in data_manager.reflection_table.keys()

  #test for successful escape from normalised intensity calculation
  n = (data_manager.reflection_table['Esq'] == 0.0).count(True)
  assert n != len(data_manager.reflection_table)

def test_KBDataManager():
  (test_reflections, test_experiments, params) = generate_test_input()
  data_manager = KB_Data_Manager(test_reflections, test_experiments, params)

def test_AimlessDataManager():
  (test_reflections, test_experiments, params) = generate_test_input()
  data_manager = AimlessDataManager(test_reflections, test_experiments, params)


def test_targeted_data_manager():
  (test_reflections, test_experiments, params) = generate_test_input()
  targeted_reflections = copy.deepcopy(test_reflections)
  targeted_reflections['inverse_scale_factor'] = flex.double([1.0, 1.0, 1.0])
  targeted_dm = TargetedDataManager(test_reflections, test_experiments,
    targeted_reflections, params)


def test_multi_data_manager():
  (test_reflections, test_experiments, params) = generate_test_input()
  params.scaling_options.__inject__('multi_mode', True)
  data_manager = MultiCrystalDataManager([test_reflections, test_reflections],
    [test_experiments, test_experiments], params)

def test_apm():
  (test_reflections, test_experiments, params) = generate_test_input()
  #params.parameterisation.absorption_term=False
  params.scaling_options.__inject__('multi_mode', False)
  data_manager = AimlessDataManager(test_reflections, test_experiments, params)
  apm = active_parameter_manager(data_manager, ['g_decay', 'g_scale'])
  assert 'g_decay' in apm.active_parameterisation
  assert 'g_scale' in apm.active_parameterisation
  assert 'g_absorption' not in apm.active_parameterisation
  assert apm.n_active_params == data_manager.g_scale.n_params + data_manager.g_decay.n_params
  assert apm.constant_g_values is not None
