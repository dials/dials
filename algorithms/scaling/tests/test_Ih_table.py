from target_Ih import SingleIhTable, JointIhTable
from dials.array_family import flex
from dials.util.options import OptionParser
from data_manager_functions import ScalingDataManager
from libtbx import phil
from dxtbx.model.experiment_list import ExperimentList
from dxtbx.model import Crystal, Scan, Beam, Goniometer
from scitbx import sparse

def generate_refl_1():
  reflections = flex.reflection_table()
  reflections['intensity.prf.value'] = flex.double(
    [100.0, 100.0, 80.0, 60.0, 30.0, 40.0, 60.0])
  reflections['intensity.prf.variance'] = flex.double(
    [90.0, 100.0, 90.0, 60.0, 30.0, 50.0, 50.0])
  reflections['miller_index'] = flex.miller_index(
    [(1, 0, 0), (0, 0, 1), (-1, 0, 0), (0, 2, 0), (0, 4, 0), (0, 0, -2), (0, 0, 2)])
  reflections['d'] = flex.double([5.0, 5.0, 5.0, 5.0, 2.5, 2.5, 2.5])
  reflections['lp'] = flex.double([1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0])
  reflections['dqe'] = flex.double([1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0])
  #reflections['xyzobs.px.value'] = flex.vec3_double([(0.0,0.0,0.0), (0.0,0.0,0.0),(0.0,0.0,0.0),
  #  (0.0,0.0,0.0), (0.0,0.0,0.0),(0.0,0.0,0.0),(0.0,0.0,0.0)])
  #reflections['s1'] = flex.vec3_double([(0.0,0.1,1.0), (0.0,0.1,1.0),(0.0,0.1,1.0),
  #  (0.0,0.1,1.0), (0.0,0.1,1.0), (0.0,0.1,1.0), (0.0,0.1,1.0)])
  reflections.set_flags(flex.bool([True, True, True, True, True, True, True]),
    reflections.flags.integrated)
  return reflections

def generate_refl_2():
  reflections = flex.reflection_table()
  reflections['intensity.prf.value'] = flex.double([60.0, 30.0])
  reflections['intensity.prf.variance'] = flex.double([60.0, 30.0])
  reflections['miller_index'] = flex.miller_index([(1, 0, 0), (0, 4, 0)])
  reflections['d'] = flex.double([5.0, 2.5])
  reflections['lp'] = flex.double([1.0, 1.0])
  reflections['dqe'] = flex.double([1.0, 1.0])
  reflections.set_flags(flex.bool([True, True]), reflections.flags.integrated)
  return reflections

def generate_test_experiments():
  experiments = ExperimentList()
  exp_dict = {"__id__" : "crystal", "real_space_a": [5.0, 0.0, 0.0],
              "real_space_b": [0.0, 10.0, 0.0], "real_space_c": [0.0, 0.0, 5.0],
              "space_group_hall_symbol": " C 2y"}
  experiments.crystal = Crystal.from_dict(exp_dict)
  experiments.scan = Scan(image_range=[0, 90], oscillation=[0.0, 1.0])
  experiments.beam = Beam(s0=(0.0, 0.0, 1.01))
  experiments.goniometer = Goniometer((1.0, 0.0, 0.0))
  return experiments

def generate_test_params():
  phil_scope = phil.parse('''
      include scope dials.algorithms.scaling.scaling_options.phil_scope
  ''', process_includes=True)
  optionparser = OptionParser(phil=phil_scope, check_format=False)
  parameters, _ = optionparser.parse_args(args=None, quick_parse=True,
    show_diff_phil=False)
  return parameters

def generate_test_datamanager():
  reflections_1 = generate_refl_1()
  experiments = generate_test_experiments()
  params = generate_test_params()
  dm = ScalingDataManager(reflections_1, experiments, params)
  return dm

def generate_second_test_datamanager():
  reflections_2 = generate_refl_2()
  experiments = generate_test_experiments()
  params = generate_test_params()
  dm = ScalingDataManager(reflections_2, experiments, params)
  return dm

def generate_single_test_input():
  dm = generate_test_datamanager()
  weights=flex.double([1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0])
  return dm.reflection_table, weights

def generate_joint_test_input():
  data_manager = generate_test_datamanager()
  weights = flex.double([1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0])
  data_manager_2 = generate_second_test_datamanager()
  weights_2 = flex.double([1.0, 1.0])
  data_manager._Ih_table = SingleIhTable(data_manager.reflection_table, weights)
  data_manager_2._Ih_table = SingleIhTable(data_manager_2.reflection_table, weights_2)
  return data_manager, data_manager_2

def test_Ih_table((reflection_table, weights)=generate_single_test_input()):
  Ih_table = SingleIhTable(reflection_table, weights)

  #upon initialisation, Ih_table should set unity scale factors, sort by miller index
  #and calculate the Ih_values. It should also create the h_index arrays and h_index_matrix
  assert (Ih_table.inverse_scale_factors == 1.0).count(True) == Ih_table.size
  assert list(Ih_table.asu_miller_index) == list(flex.miller_index(
    [(0, 0, 1), (0, 0, 2), (0, 0, 2), (0, 2, 0), (0, 4, 0), (1, 0, 0), (1, 0, 0)]))
  assert list(Ih_table.h_index_counter_array) == list(flex.int([1, 2, 1, 1, 2]))
  assert list(Ih_table.h_index_cumulative_array) == list(flex.int([0, 1, 3, 4, 5, 7]))
  assert list(Ih_table.n_h) == list(flex.double([1, 2, 2, 1, 1, 2, 2]))
  assert list(Ih_table.Ih_values) == list(flex.double(
    [100.0, 50.0, 50.0, 60.0, 30.0, 90.0, 90.0]))

  assert Ih_table.h_index_matrix[0, 0] == 1
  assert Ih_table.h_index_matrix[1, 1] == 1
  assert Ih_table.h_index_matrix[2, 1] == 1
  assert Ih_table.h_index_matrix[3, 2] == 1
  assert Ih_table.h_index_matrix[4, 3] == 1
  assert Ih_table.h_index_matrix[5, 4] == 1
  assert Ih_table.h_index_matrix[6, 4] == 1
  assert Ih_table.h_index_matrix.non_zeroes == 7
  assert Ih_table.h_index_matrix.n_cols == 5
  assert Ih_table.h_index_matrix.n_rows == 7


def test_joint_Ih_table((dm1, dm2)=generate_joint_test_input()):
  Ih_table = JointIhTable([dm1, dm2])

  #test that the values have been sorted/combined correctly
  assert list(Ih_table.asu_miller_index) == list(flex.miller_index([(0 ,0, 1),
    (0, 0, 2), (0, 0, 2), (0, 2, 0), (0, 4, 0), (0, 4, 0), (1, 0, 0), (1, 0, 0),
    (1, 0, 0)]))
  assert list(Ih_table.Ih_values) == list(flex.double(
    [100.0, 50.0, 50.0, 60.0, 30.0, 30.0, 80.0, 80.0, 80.0]))
  assert list(Ih_table.h_index_counter_array) == list(flex.int([1, 2, 1, 2, 3]))
  assert list(Ih_table.h_index_cumulative_array) == list(flex.int([0, 1, 3, 4, 6, 9]))
  assert list(Ih_table._h_idx_count_list[0]) == list(flex.int([1, 2, 1, 1, 2]))
  assert list(Ih_table._h_idx_count_list[1]) == list(flex.int([0, 0, 0, 1, 1]))
  assert list(Ih_table._h_idx_cumulative_list[0]) == list(flex.int([0, 1, 3, 4, 5, 7]))
  assert list(Ih_table._h_idx_cumulative_list[1]) == list(flex.int([0, 0, 0, 0, 1, 2]))

  #now test h_expand matrices
  assert Ih_table.h_index_expand_list[0][0, 0] == 1
  assert Ih_table.h_index_expand_list[0][1, 1] == 1
  assert Ih_table.h_index_expand_list[0][2, 2] == 1
  assert Ih_table.h_index_expand_list[0][3, 3] == 1
  assert Ih_table.h_index_expand_list[0][4, 4] == 1
  assert Ih_table.h_index_expand_list[0][5, 6] == 1
  assert Ih_table.h_index_expand_list[0][6, 7] == 1
  assert Ih_table.h_index_expand_list[0].non_zeroes == 7
  assert Ih_table.h_index_expand_list[0].n_cols == 9
  assert Ih_table.h_index_expand_list[0].n_rows == 7

  assert Ih_table.h_index_expand_list[1][0, 5] == 1
  assert Ih_table.h_index_expand_list[1][1, 8] == 1
  assert Ih_table.h_index_expand_list[1].non_zeroes == 2
  assert Ih_table.h_index_expand_list[1].n_cols == 9
  assert Ih_table.h_index_expand_list[1].n_rows == 2
