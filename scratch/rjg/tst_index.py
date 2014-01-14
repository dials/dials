from __future__ import division
try:
  import scipy.linalg # import dependency
except ImportError, e:
  pass

import os

import libtbx.load_env
from libtbx import easy_run
from libtbx.test_utils import approx_equal
from libtbx.test_utils import open_tmp_directory
from scitbx import matrix
from cctbx.array_family import flex # import dependency
from cctbx import uctbx
from dials.model.serialize import load

have_dials_regression = libtbx.env.has_module("dials_regression")
if have_dials_regression:
  dials_regression = libtbx.env.find_in_repositories(
    relative_path="dials_regression",
    test=os.path.isdir)

# check optional dependencies of multi-lattice search
try:
  import scipy
  have_scipy = True
except ImportError:
  have_scipy = False

try:
  import sklearn
  have_sklearn = True
except ImportError:
  have_sklearn = False

try:
  import hcluster
  have_hcluster = True
except ImportError:
  have_hcluster = False

try:
  import networkx
  have_networkx = True
except ImportError:
  have_networkx = False


def check_external_dependencies(dependencies):
  missing = []
  for dependency in dependencies:
    if not globals().get('have_%s' %dependency, False):
      missing.append(dependency)
  return missing


class run_one_indexing(object):
  def __init__(self,
               pickle_path,
               sweep_path,
               extra_args,
               expected_unit_cell,
               expected_rmsds,
               expected_hall_symbol,
               n_expected_lattices=1,
               relative_length_tolerance=0.005,
               absolute_angle_tolerance=0.5):

    py_script_path = libtbx.env.find_in_repositories(
      relative_path="dials/scratch/rjg/index_3D_FFT_simple.py",
      test=os.path.isfile)

    args = ["dials.python", py_script_path, pickle_path, sweep_path] + extra_args

    cwd = os.path.abspath(os.curdir)
    tmp_dir = open_tmp_directory(suffix="test_3DFFT_indexing")
    os.chdir(tmp_dir)
    command = " ".join(args)
    result = easy_run.fully_buffered(command=command).raise_if_errors()
    os.chdir(cwd)
    for i in range(n_expected_lattices):
      suffix = ""
      if n_expected_lattices > 1:
        suffix = "_%i" %(i+1)
      assert os.path.exists(os.path.join(tmp_dir, "indexed%s.pickle" %suffix))
      assert os.path.exists(os.path.join(tmp_dir, "crystal%s.json" %suffix))
      assert os.path.exists(os.path.join(tmp_dir, "sweep%s.json" %suffix))
      self.reflections = load.reflections(os.path.join(tmp_dir, "indexed%s.pickle" %suffix))
      self.crystal_model = load.crystal(os.path.join(tmp_dir, "crystal%s.json" %suffix))
      self.sweep = load.sweep(os.path.join(tmp_dir, "sweep%s.json" %suffix))
      assert self.crystal_model.get_unit_cell().is_similar_to(
        expected_unit_cell,
        relative_length_tolerance=relative_length_tolerance,
        absolute_angle_tolerance=absolute_angle_tolerance)
      sg = self.crystal_model.get_space_group()
      assert sg.type().hall_symbol() == expected_hall_symbol
      mi = self.reflections.miller_index()
      assert (mi != (0,0,0)).count(False) == 0
      self.reflections = self.reflections.select(mi != (0,0,0))
      self.rmsds = self.get_rmsds_obs_pred(
        self.reflections, self.sweep, self.crystal_model)
      for actual, expected in zip(self.rmsds, expected_rmsds):
        assert actual <= expected

  def get_rmsds_obs_pred(self, observations, sweep, crystal_model):
    from dials.algorithms.spot_prediction import ray_intersection
    reflections = ray_intersection(sweep.get_detector(), observations)
    from dials.scratch.rjg.index_3D_FFT_simple import master_params
    from dials.algorithms.refinement import RefinerFactory
    refine = RefinerFactory.from_parameters_data_models(
      master_params, reflections, sweep=sweep, crystal=crystal_model,
      verbosity=0)
    return refine.rmsds()

def exercise_1():
  # thaumatin
  data_dir = os.path.join(dials_regression, "indexing_test_data", "i04_weak_data")
  pickle_path = os.path.join(data_dir, "full.pickle")
  sweep_path = os.path.join(data_dir, "sweep_orig.json")
  extra_args = ["multiple_lattice_search=False", # use older non-clustering version
                "reflections_per_degree=5",
                "n_macro_cycles=2"]
  expected_unit_cell = uctbx.unit_cell(
    (58, 58, 150, 90, 90, 90))
  expected_rmsds = (0.06, 0.05, 0.0005)
  expected_hall_symbol = ' P 1'

  result = run_one_indexing(pickle_path, sweep_path, extra_args, expected_unit_cell,
                     expected_rmsds, expected_hall_symbol)

def exercise_2():
  missing = check_external_dependencies(['scipy', 'sklearn', 'networkx'])
  if len(missing):
    print ("Skipping exercise_2: missing dependencies" +
           " %s" * len(missing)) %(tuple(missing))
    return
  # thaumatin
  data_dir = os.path.join(dials_regression, "indexing_test_data", "i04_weak_data")
  pickle_path = os.path.join(data_dir, "full.pickle")
  sweep_path = os.path.join(data_dir, "sweep_orig.json")
  extra_args = ["multiple_lattice_search=True",
                "reflections_per_degree=5",
                "n_macro_cycles=2",
                "bin_size_fraction=0.25",
                "d_min=4"]
  expected_unit_cell = uctbx.unit_cell(
    (58, 58, 150, 90, 90, 90))
  expected_rmsds = (0.06, 0.05, 0.0005)
  expected_hall_symbol = ' P 1'

  result = run_one_indexing(pickle_path, sweep_path, extra_args, expected_unit_cell,
                            expected_rmsds, expected_hall_symbol)

def exercise_3():
  missing = check_external_dependencies(['scipy', 'sklearn', 'networkx'])
  if len(missing):
    print ("Skipping exercise_2: missing dependencies" +
           " %s" * len(missing)) %(tuple(missing))
    return
  # thaumatin
  data_dir = os.path.join(dials_regression, "indexing_test_data", "i04_weak_data")
  pickle_path = os.path.join(data_dir, "full.pickle")
  sweep_path = os.path.join(data_dir, "sweep_orig.json")
  extra_args = ["multiple_lattice_search=True",
                "reflections_per_degree=5",
                "n_macro_cycles=2",
                "d_min=4"]
  expected_unit_cell = uctbx.unit_cell(
    (58, 58, 150, 90, 90, 90))
  expected_rmsds = (0.06, 0.05, 0.0005)

  # now enforce symmetry
  extra_args.append("space_group=P4")
  expected_hall_symbol = ' P 4'

  result = run_one_indexing(pickle_path, sweep_path, extra_args, expected_unit_cell,
                            expected_rmsds, expected_hall_symbol)

  a, b, c = result.crystal_model.get_real_space_vectors()
  assert approx_equal(a.length(), b.length())
  assert c.length() > b.length()
  assert approx_equal(a.angle(b, deg=True), 90)
  assert approx_equal(b.angle(c, deg=True), 90)
  assert approx_equal(c.angle(a, deg=True), 90)

def exercise_4():
  missing = check_external_dependencies(['scipy', 'sklearn', 'networkx'])
  if len(missing):
    print ("Skipping exercise_2: missing dependencies" +
           " %s" * len(missing)) %(tuple(missing))
    return
  # trypsin
  data_dir = os.path.join(dials_regression, "indexing_test_data", "trypsin")
  pickle_path = os.path.join(data_dir, "P1_X6_1.pickle")
  sweep_path = os.path.join(data_dir, "sweep_P1_X6_1.json")
  extra_args = ["multiple_lattice_search=True",
                "reflections_per_degree=5",
                "n_macro_cycles=2",
                "d_min=4",
                "scan_range=0,50",
                "scan_range=450,500",
                "scan_range=850,900"]
  expected_unit_cell = uctbx.unit_cell(
    (54.3, 58.3, 66.5, 90, 90, 90))
  expected_rmsds = (0.08, 0.06, 0.002)
  expected_hall_symbol = ' P 1'

  result = run_one_indexing(pickle_path, sweep_path, extra_args, expected_unit_cell,
                            expected_rmsds, expected_hall_symbol)

def exercise_5():
  missing = check_external_dependencies(['scipy', 'sklearn', 'networkx'])
  if len(missing):
    print ("Skipping exercise_2: missing dependencies" +
           " %s" * len(missing)) %(tuple(missing))
    return
  # synthetic trypsin multi-lattice dataset (2 lattices)
  data_dir = os.path.join(dials_regression, "indexing_test_data", "trypsin")
  pickle_path = os.path.join(data_dir, "P1_X6_1_2.pickle")
  sweep_path = os.path.join(data_dir, "sweep_P1_X6_1_2.json")
  extra_args = ["multiple_lattice_search=True",
                "reflections_per_degree=5",
                "n_macro_cycles=2",
                "d_min=4",
                "scan_range=0,50",
                "scan_range=450,500",
                "scan_range=850,900"]
  expected_unit_cell = uctbx.unit_cell(
    (54.3, 58.3, 66.5, 90, 90, 90))
  expected_rmsds = (0.17, 0.20, 0.004)
  expected_hall_symbol = ' P 1'
  n_expected_lattices = 2

  result = run_one_indexing(pickle_path, sweep_path, extra_args, expected_unit_cell,
                            expected_rmsds, expected_hall_symbol,
                            n_expected_lattices=n_expected_lattices,
                            relative_length_tolerance=0.02,
                            absolute_angle_tolerance=1)

def exercise_6():
  missing = check_external_dependencies(['scipy', 'sklearn', 'networkx'])
  if len(missing):
    print ("Skipping exercise_2: missing dependencies" +
           " %s" * len(missing)) %(tuple(missing))
    return
  # synthetic trypsin multi-lattice dataset (3 lattices)
  data_dir = os.path.join(dials_regression, "indexing_test_data", "trypsin")
  pickle_path = os.path.join(data_dir, "P1_X6_1_2_3.pickle")
  sweep_path = os.path.join(data_dir, "sweep_P1_X6_1_2_3.json")
  extra_args = ["multiple_lattice_search=True",
                "reflections_per_degree=5",
                "n_macro_cycles=2",
                "d_min=4",
                "max_cell=70", #XXX eventually this should not be needed
                ]
  expected_unit_cell = uctbx.unit_cell(
    (54.3, 58.3, 66.5, 90, 90, 90))
  expected_rmsds = (0.24, 0.30, 0.006)
  expected_hall_symbol = ' P 1'
  n_expected_lattices = 3

  result = run_one_indexing(pickle_path, sweep_path, extra_args, expected_unit_cell,
                            expected_rmsds, expected_hall_symbol,
                            n_expected_lattices=n_expected_lattices,
                            relative_length_tolerance=0.01,
                            absolute_angle_tolerance=1)

def exercise_7():
  missing = check_external_dependencies(['scipy', 'sklearn', 'networkx'])
  if len(missing):
    print ("Skipping exercise_2: missing dependencies" +
           " %s" * len(missing)) %(tuple(missing))
    return
  # synthetic trypsin multi-lattice dataset (4 lattices)
  data_dir = os.path.join(dials_regression, "indexing_test_data", "trypsin")
  pickle_path = os.path.join(data_dir, "P1_X6_1_2_3_4.pickle")
  sweep_path = os.path.join(data_dir, "sweep_P1_X6_1_2_3_4.json")
  extra_args = ["multiple_lattice_search=True",
                "reflections_per_degree=5",
                "n_macro_cycles=2",
                "d_min=4",
                "max_cell=70", #XXX eventually this should not be needed
                ]
  expected_unit_cell = uctbx.unit_cell(
    (54.3, 58.3, 66.5, 90, 90, 90))
  expected_rmsds = (0.24, 0.30, 0.006)
  expected_hall_symbol = ' P 1'
  n_expected_lattices = 4

  result = run_one_indexing(pickle_path, sweep_path, extra_args, expected_unit_cell,
                            expected_rmsds, expected_hall_symbol,
                            n_expected_lattices=n_expected_lattices,
                            relative_length_tolerance=0.01,
                            absolute_angle_tolerance=1)

def exercise_8():
  # synthetic trypsin multi-lattice dataset (4 lattices)
  data_dir = os.path.join(dials_regression, "indexing_test_data", "trypsin")
  pickle_path = os.path.join(data_dir, "P1_X6_1_2_3_4.pickle")
  sweep_path = os.path.join(data_dir, "sweep_P1_X6_1_2_3_4.json")
  extra_args = ["method=real_space_grid_search",
                "use_all_reflections=True",
                "n_macro_cycles=5",
                "d_min=4",
                "known_symmetry.unit_cell=54.3,58.3,66.5,90,90,90",
                "known_symmetry.space_group=P212121",
                "scan_range=0,10",
                "engine=LevMarIterations",
                "beam.fix=all",
                "detector.fix=orientation",
                "maximum_spot_error=3",
                #"max_cell=70", #XXX eventually this should not be needed
                ]
  expected_unit_cell = uctbx.unit_cell((54.3, 58.3, 66.5, 90, 90, 90))
  expected_rmsds = (0.28, 0.30, 0.006)
  expected_hall_symbol = ' P 2ac 2ab'
  n_expected_lattices = 1

  result = run_one_indexing(pickle_path, sweep_path, extra_args, expected_unit_cell,
                            expected_rmsds, expected_hall_symbol,
                            n_expected_lattices=n_expected_lattices,
                            relative_length_tolerance=0.02,
                            absolute_angle_tolerance=1)


def run(args):
  if not libtbx.env.has_module("dials_regression"):
    print "Skipping exercise_index_3D_FFT_simple: dials_regression not present"
    return

  exercises = (exercise_1, exercise_2, exercise_3, exercise_4, exercise_5,
               exercise_6, exercise_7, exercise_8)
  if len(args):
    args = [int(arg) for arg in args]
    for arg in args: assert arg > 0
    exercises = [exercises[arg-1] for arg in args]

  for exercise in exercises:
    exercise()

if __name__ == '__main__':
  import sys
  from libtbx.utils import show_times_at_exit
  show_times_at_exit()
  run(sys.argv[1:])
