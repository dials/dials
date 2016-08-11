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
from cctbx.array_family import flex # import dependency
from cctbx import uctbx
from dials.model.serialize import load
from dxtbx.serialize import load as dxtbx_load

have_dials_regression = libtbx.env.has_module("dials_regression")
if have_dials_regression:
  dials_regression = libtbx.env.find_in_repositories(
    relative_path="dials_regression",
    test=os.path.isdir)

# check optional dependencies of multi-lattice search
try:
  import scipy # implicit import
  have_scipy = True
except ImportError:
  have_scipy = False

try:
  import sklearn # implicit import
  have_sklearn = True
except ImportError:
  have_sklearn = False

try:
  import networkx # implicit import
  have_networkx = True
except ImportError:
  have_networkx = False


def check_external_dependencies(dependencies):
  missing = []
  for dependency in dependencies:
    if not globals().get('have_%s' %dependency, False):
      missing.append(dependency)
  return missing

def unit_cells_are_similar(uc1, uc2, relative_length_tolerance=0.01,
                           absolute_angle_tolerance=1):
  # see also uctbx.cpp unit_cell::is_similar_to()
  l1 = uc1.parameters()
  l2 = uc2.parameters()
  for i in range(3):
    if abs(min(l1[i], l2[i]) / max(l1[i], l2[i]) - 1) > relative_length_tolerance:
      return False
  for i in range(3, 6):
    if abs(l1[i]-l2[i]) > absolute_angle_tolerance:
      if abs(l1[i]-(180-l2[i])) > absolute_angle_tolerance:
        return False
      #else:
        #print uc1, uc2
  return True

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

    args = ["dials.index", pickle_path, sweep_path] + extra_args

    cwd = os.path.abspath(os.curdir)
    tmp_dir = open_tmp_directory(suffix="test_dials_index")
    os.chdir(tmp_dir)
    command = " ".join(args)
    print command
    result = easy_run.fully_buffered(command=command).raise_if_errors()
    os.chdir(cwd)
    assert os.path.exists(os.path.join(tmp_dir, "experiments.json"))
    experiments_list = dxtbx_load.experiment_list(
      os.path.join(tmp_dir, "experiments.json"), check_format=False)
    assert len(experiments_list.crystals()) == n_expected_lattices, (
      len(experiments_list.crystals()), n_expected_lattices)
    assert os.path.exists(os.path.join(tmp_dir, "indexed.pickle"))
    from libtbx.utils import time_log
    unpickling_timer = time_log("unpickling")
    self.calc_rmsds_timer = time_log("calc_rmsds")
    unpickling_timer.start()
    self.indexed_reflections = load.reflections(os.path.join(tmp_dir, "indexed.pickle"))
    unpickling_timer.stop()
    for i in range(len(experiments_list)):
      experiment = experiments_list[i]
      self.crystal_model = experiment.crystal
      #assert self.crystal_model.get_unit_cell().is_similar_to(
        #expected_unit_cell,
        #relative_length_tolerance=relative_length_tolerance,
        #absolute_angle_tolerance=absolute_angle_tolerance), (
          #self.crystal_model.get_unit_cell().parameters(),
          #expected_unit_cell.parameters())
      assert unit_cells_are_similar(
        self.crystal_model.get_unit_cell(),expected_unit_cell,
        relative_length_tolerance=relative_length_tolerance,
        absolute_angle_tolerance=absolute_angle_tolerance), (
          self.crystal_model.get_unit_cell().parameters(),
          expected_unit_cell.parameters())
      sg = self.crystal_model.get_space_group()
      assert sg.type().hall_symbol() == expected_hall_symbol, (
        sg.type().hall_symbol(), expected_hall_symbol)
      reflections = self.indexed_reflections.select(
        self.indexed_reflections['id'] == i)
      mi = reflections['miller_index']
      assert (mi != (0,0,0)).count(False) == 0
      reflections = reflections.select(mi != (0,0,0))
      self.rmsds = self.get_rmsds_obs_pred(reflections, experiment)
      for actual, expected in zip(self.rmsds, expected_rmsds):
        assert actual <= expected, "%s %s" %(self.rmsds, expected_rmsds)
    if 0:
      print self.calc_rmsds_timer.legend
      print unpickling_timer.report()
      print self.calc_rmsds_timer.report()

  def get_rmsds_obs_pred(self, observations, experiment):
    reflections = observations.select(observations.get_flags(
      observations.flags.used_in_refinement))
    assert len(reflections) > 0
    obs_x, obs_y, obs_z = reflections['xyzobs.mm.value'].parts()
    calc_x, calc_y, calc_z = reflections['xyzcal.mm'].parts()
    rmsd_x = flex.mean(flex.pow2(obs_x-calc_x))**0.5
    rmsd_y = flex.mean(flex.pow2(obs_y-calc_y))**0.5
    rmsd_z = flex.mean(flex.pow2(obs_z-calc_z))**0.5
    return (rmsd_x, rmsd_y, rmsd_z)

def exercise_1():
  # thaumatin
  data_dir = os.path.join(dials_regression, "indexing_test_data", "i04_weak_data")
  pickle_path = os.path.join(data_dir, "full.pickle")
  sweep_path = os.path.join(data_dir, "datablock_orig.json")
  extra_args = ["bin_size_fraction=0.25",
                "scan_range=1,20",
                "scan_range=250,270",
                "scan_range=520,540"]
  expected_unit_cell = uctbx.unit_cell(
    (58, 58, 150, 90, 90, 90))
  expected_rmsds = (0.05, 0.04, 0.0005)
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
  sweep_path = os.path.join(data_dir, "datablock_orig.json")
  extra_args = ["cluster_analysis_search=True",
                "n_macro_cycles=3",
                "bin_size_fraction=0.25",
                "reciprocal_space_grid.d_min=4"]
  expected_unit_cell = uctbx.unit_cell(
    (58, 58, 150, 90, 90, 90))
  expected_rmsds = (0.05, 0.04, 0.0004)
  expected_hall_symbol = ' P 1'

  result = run_one_indexing(pickle_path, sweep_path, extra_args, expected_unit_cell,
                            expected_rmsds, expected_hall_symbol)

def exercise_3():
  missing = check_external_dependencies(['scipy', 'sklearn', 'networkx'])
  if len(missing):
    print ("Skipping exercise_3: missing dependencies" +
           " %s" * len(missing)) %(tuple(missing))
    return
  # thaumatin
  data_dir = os.path.join(dials_regression, "indexing_test_data", "i04_weak_data")
  pickle_path = os.path.join(data_dir, "full.pickle")
  sweep_path = os.path.join(data_dir, "datablock_orig.json")
  extra_args = ["cluster_analysis_search=True",
                "n_macro_cycles=3",
                "bin_size_fraction=0.25",
                "reciprocal_space_grid.d_min=4"]
  expected_unit_cell = uctbx.unit_cell(
    (58, 58, 150, 90, 90, 90))
  expected_rmsds = (0.05, 0.041, 0.0004)

  # now enforce symmetry
  extra_args.append("known_symmetry.space_group=P4")
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
    print ("Skipping exercise_4: missing dependencies" +
           " %s" * len(missing)) %(tuple(missing))
    return
  # trypsin
  data_dir = os.path.join(dials_regression, "indexing_test_data", "trypsin")
  pickle_path = os.path.join(data_dir, "P1_X6_1.pickle")
  sweep_path = os.path.join(data_dir, "datablock_P1_X6_1.json")
  extra_args = ["cluster_analysis_search=True",
                "n_macro_cycles=3",
                "reciprocal_space_grid.d_min=4",
                "filter_overlaps=False", # P1_X6_1.pickle does not contain bbox!
                "scan_range=0,50",
                "scan_range=450,500",
                "scan_range=850,900"]
  expected_unit_cell = uctbx.unit_cell(
    (54.3, 58.3, 66.5, 90, 90, 90))
  expected_rmsds = (0.061, 0.06, 0.00042)
  expected_hall_symbol = ' P 1'

  result = run_one_indexing(pickle_path, sweep_path, extra_args, expected_unit_cell,
                            expected_rmsds, expected_hall_symbol)

def exercise_5():
  missing = check_external_dependencies(['scipy', 'sklearn', 'networkx'])
  if len(missing):
    print ("Skipping exercise_5: missing dependencies" +
           " %s" * len(missing)) %(tuple(missing))
    return
  # synthetic trypsin multi-lattice dataset (2 lattices)
  data_dir = os.path.join(dials_regression, "indexing_test_data", "trypsin")
  pickle_path = os.path.join(data_dir, "P1_X6_1_2.pickle")
  sweep_path = os.path.join(data_dir, "datablock_P1_X6_1_2.json")
  extra_args = ["cluster_analysis_search=True",
                "reflections_per_degree=10",
                "n_macro_cycles=2",
                "reciprocal_space_grid.d_min=4",
                "max_cell=70",
                "scan_range=0,50",
                "scan_range=450,500",
                "scan_range=850,900"]
  expected_unit_cell = uctbx.unit_cell(
    (54.3, 58.3, 66.5, 90, 90, 90))
  expected_rmsds = (0.1, 0.1, 0.005)
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
    print ("Skipping exercise_6: missing dependencies" +
           " %s" * len(missing)) %(tuple(missing))
    return
  # synthetic trypsin multi-lattice dataset (3 lattices)
  data_dir = os.path.join(dials_regression, "indexing_test_data", "trypsin")
  pickle_path = os.path.join(data_dir, "P1_X6_1_2_3.pickle")
  sweep_path = os.path.join(data_dir, "datablock_P1_X6_1_2_3.json")
  extra_args = ["cluster_analysis_search=True",
                "reflections_per_degree=10",
                "n_macro_cycles=2",
                "reciprocal_space_grid.d_min=4",
                "max_cell=70",
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
    print ("Skipping exercise_7: missing dependencies" +
           " %s" * len(missing)) %(tuple(missing))
    return
  # synthetic trypsin multi-lattice dataset (4 lattices)
  data_dir = os.path.join(dials_regression, "indexing_test_data", "trypsin")
  pickle_path = os.path.join(data_dir, "P1_X6_1_2_3_4.pickle")
  sweep_path = os.path.join(data_dir, "datablock_P1_X6_1_2_3_4.json")
  extra_args = ["cluster_analysis_search=True",
                "reflections_per_degree=10",
                "n_macro_cycles=2",
                "reciprocal_space_grid.d_min=4",
                "max_cell=70",
                "max_lattices=4"
                ]
  expected_unit_cell = uctbx.unit_cell(
    (54.3, 58.3, 66.5, 90, 90, 90))
  expected_rmsds = (0.24, 0.23, 0.006)
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
  sweep_path = os.path.join(data_dir, "datablock_P1_X6_1_2_3_4.json")
  extra_args = ["indexing.method=real_space_grid_search",
                "reflections_per_degree=10",
                "n_macro_cycles=5",
                "known_symmetry.unit_cell=54.3,58.3,66.5,90,90,90",
                "known_symmetry.space_group=P212121",
                "scan_range=0,10",
                "beam.fix=all",
                "detector.fix=all",
                "max_cell=70",
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

def exercise_9():
  # thaumatin
  data_dir = os.path.join(dials_regression, "indexing_test_data", "i04_weak_data")
  pickle_path = os.path.join(data_dir, "full.pickle")
  sweep_path = os.path.join(data_dir, "datablock_orig.json")
  extra_args = ["n_macro_cycles=2",
                "indexing.method=fft1d",
                "bin_size_fraction=0.25",
                "scan_range=1,20",
                "scan_range=250,270",
                "scan_range=520,540"]
  expected_unit_cell = uctbx.unit_cell(
    (58, 58, 150, 90, 90, 90))
  expected_rmsds = (0.06, 0.05, 0.0005)
  expected_hall_symbol = ' P 1'

  result = run_one_indexing(pickle_path, sweep_path, extra_args, expected_unit_cell,
                     expected_rmsds, expected_hall_symbol)

def exercise_10():
  # synthetic trypsin multi-lattice dataset (3 lattices)
  data_dir = os.path.join(dials_regression, "indexing_test_data", "trypsin")
  pickle_path = os.path.join(data_dir, "P1_X6_1_2_3.pickle")
  sweep_path = os.path.join(data_dir, "datablock_P1_X6_1_2_3.json")
  extra_args = ["indexing.method=real_space_grid_search",
                "d_min_start=3",
                "n_macro_cycles=3",
                "known_symmetry.unit_cell=54.3,58.3,66.5,90,90,90",
                "known_symmetry.space_group=P212121",
                "scan_range=0,10",
                "beam.fix=all",
                "detector.fix=all",
                "max_lattices=3",
                "index_assignment.method=local",
                "nearest_neighbours=50",
                ]

  expected_unit_cell = uctbx.unit_cell((54.3, 58.3, 66.5, 90, 90, 90))
  expected_rmsds = (0.33, 0.40, 0.0022)
  expected_hall_symbol = ' P 2ac 2ab'
  n_expected_lattices = 3

  result = run_one_indexing(pickle_path, sweep_path, extra_args, expected_unit_cell,
                            expected_rmsds, expected_hall_symbol,
                            n_expected_lattices=n_expected_lattices,
                            relative_length_tolerance=0.02,
                            absolute_angle_tolerance=1)

def exercise_11():
  return # disable test until a better image and/or parameters are available
#  image_path = os.path.join(dials_regression, "spotfinding_test_data",
#                            "idx-s00-20131106040304531.cbf")
#  cwd = os.path.abspath(os.curdir)
#  tmp_dir = os.path.abspath(open_tmp_directory(suffix="test_dials_index"))
#  os.chdir(tmp_dir)
#  print tmp_dir
#
#  args = ["dials.import", image_path,
#          "output.datablock=datablock.json"]
#  command = " ".join(args)
#  #print command
#  result = easy_run.fully_buffered(command=command).raise_if_errors()
#
#  datablock_json = os.path.join(tmp_dir, "datablock.json")
#
#  args = ["dials.find_spots",
#          datablock_json,
#          "threshold.xds.sigma_strong=7",
#          "min_spot_size=6",
#          ]
#
#  command = " ".join(args)
#  #print command
#  result = easy_run.fully_buffered(command=command).raise_if_errors()
#  pickle_path = os.path.join(tmp_dir, "strong.pickle")
#  assert os.path.exists(pickle_path)
#
#  extra_args = ["indexing.method=real_space_grid_search",
#                "n_macro_cycles=3",
#                "known_symmetry.unit_cell=78,78,39,90,90,90",
#                "known_symmetry.space_group=P43212",
#                "beam.fix=all",
#                "detector.fix=all",
#                "hkl_tolerance=0.5",
#                ]
#
#  expected_unit_cell = uctbx.unit_cell((78, 78, 39, 90, 90, 90))
#  expected_rmsds = (0.31, 0.38) # XXX these rmsds really aren't great
#  expected_hall_symbol = ' P 4nw 2abw'
#  n_expected_lattices = 1
#
#  result = run_one_indexing(pickle_path, datablock_json, extra_args, expected_unit_cell,
#                            expected_rmsds, expected_hall_symbol,
#                            n_expected_lattices=n_expected_lattices,
#                            relative_length_tolerance=0.05,
#                            absolute_angle_tolerance=1)

def exercise_12():
  missing = check_external_dependencies(['scipy', 'sklearn', 'networkx'])
  if len(missing):
    print ("Skipping exercise_12: missing dependencies" +
           " %s" * len(missing)) %(tuple(missing))
    return
  # test indexing from single image of i04_weak_data
  data_dir = os.path.join(dials_regression, "indexing_test_data", "i04_weak_data")
  pickle_path = os.path.join(data_dir, "first_image.pickle")
  sweep_path = os.path.join(data_dir, "datablock_orig.json")
  extra_args = ["indexing.method=fft3d",
                "known_symmetry.space_group=P4",
                "known_symmetry.unit_cell=57.8,57.8,150,90,90,90",
                "peak_search=clean",
                "cluster_analysis_search=True",
                "min_samples=15",
                "n_macro_cycles=4",
                "reciprocal_space_grid.d_min=4"
                ]

  expected_unit_cell = uctbx.unit_cell((57.8,57.8,150,90,90,90))
  expected_rmsds = (0.06, 0.07, 0.003)
  expected_hall_symbol = ' P 4'

  result = run_one_indexing(pickle_path, sweep_path, extra_args, expected_unit_cell,
                            expected_rmsds, expected_hall_symbol)

def exercise_13():
  # test on spots derived from imosflm tutorial data:
  # http://www.ccp4.ac.uk/courses/BCA2005/tutorials/dataproc-tutorial.html
  data_dir = os.path.join(dials_regression, "indexing_test_data", "imosflm_hg_mar")
  pickle_path = os.path.join(data_dir, "strong.pickle")
  sweep_path = os.path.join(data_dir, "datablock.json")

  unit_cell = uctbx.unit_cell((58.373, 58.373, 155.939, 90, 90, 120))
  hall_symbol = '-R 3 2"'

  for uc, hall in ((unit_cell, hall_symbol), (None, hall_symbol)):
    extra_args = ["bin_size_fraction=0.25"]
    if uc is not None:
      extra_args.append("known_symmetry.unit_cell=\"%s %s %s %s %s %s\"" %unit_cell.parameters())
    if hall is not None:
      extra_args.append("known_symmetry.space_group=\"Hall: %s\"" %hall.replace('"', '\\"'))

    expected_unit_cell = unit_cell
    if hall is not None:
      expected_hall_symbol = hall
    else:
      expected_hall_symbol = ' P 1'
    expected_rmsds = (0.08, 0.11, 0.004)

    result = run_one_indexing(pickle_path, sweep_path, extra_args, expected_unit_cell,
                              expected_rmsds, expected_hall_symbol)

def exercise_14():
  from glob import glob
  data_dir = os.path.join(dials_regression, "xia2_demo_data")

  cwd = os.path.abspath(os.curdir)
  tmp_dir = os.path.abspath(open_tmp_directory())
  os.chdir(tmp_dir)
  print tmp_dir

  import shutil
  for i, image_path in enumerate(("insulin_1_001.img", "insulin_1_045.img")):
    shutil.copyfile(
      os.path.join(data_dir, image_path), "image_00%i.img" %(i+1))

  args = ["dials.import", ' '.join(glob(os.path.join(tmp_dir, "image_00*.img"))),
          "output.datablock=datablock.json", "allow_multiple_sweeps=True"]
  command = " ".join(args)
  #print command
  result = easy_run.fully_buffered(command=command).raise_if_errors()

  datablock_json = os.path.join(tmp_dir, "datablock.json")

  args = ["dials.find_spots", datablock_json]

  command = " ".join(args)
  print command
  result = easy_run.fully_buffered(command=command).raise_if_errors()
  pickle_path = os.path.join(tmp_dir, "strong.pickle")
  assert os.path.exists(pickle_path)

  expected_unit_cell = uctbx.unit_cell((78.163, 78.163, 78.163, 90.000, 90.000, 90.000))
  expected_hall_symbol = ' I 2 2 3'
  expected_rmsds = (0.05, 0.06, 0.01)

  for method in ("fft3d", "fft1d", "real_space_grid_search"):
    extra_args = []
    extra_args.append(
      "known_symmetry.unit_cell=\"%s %s %s %s %s %s\"" %expected_unit_cell.parameters())
    extra_args.append("known_symmetry.space_group=\"Hall: %s\"" %expected_hall_symbol)
    extra_args.append("indexing.method=%s" %method)
    extra_args.append("treat_single_image_as_still=False")

    result = run_one_indexing(pickle_path, datablock_json, extra_args, expected_unit_cell,
                              expected_rmsds, expected_hall_symbol)

def exercise_15():
  data_dir = os.path.join(dials_regression, "indexing_test_data", "4rotation")
  pickle_path = os.path.join(data_dir, "strong.pickle")
  sweep_path = os.path.join(data_dir, "datablock_import.json")
  extra_args = ["max_try=10", "reflections_per_degree=50",
                "known_symmetry.space_group=R3",
                "n_macro_cycles=3"]
  expected_unit_cell = uctbx.unit_cell((48.397, 48.397, 284.767, 90, 90, 120))
  expected_rmsds = (0.06, 0.08, 0.22)
  expected_hall_symbol = ' R 3'

  result = run_one_indexing(pickle_path, sweep_path, extra_args, expected_unit_cell,
                            expected_rmsds, expected_hall_symbol)
  assert len(result.indexed_reflections) > 276800, len(result.indexed_reflections)

def exercise_16():
  # test for small molecule multi-sweep indexing, 4 sweeps with different values
  # of goniometer.fixed_rotation()
  data_dir = os.path.join(dials_regression, "indexing_test_data", "multi_sweep")
  import glob
  pickle_paths = [
    glob.glob(os.path.join(data_dir, "SWEEP%i" %(i+1), "index", "*_strong.pickle"))[0]
    for i in range(4)]
  sweep_paths = [
    glob.glob(os.path.join(data_dir, "SWEEP%i" %(i+1), "index", "*_datablock_import.json"))[0]
    for i in range(4)]
  extra_args = ["known_symmetry.space_group=I4", "filter_ice=False"]
  expected_unit_cell = uctbx.unit_cell(
    (7.310, 7.310, 6.820, 90.000, 90.000, 90.000))
  expected_rmsds = (0.10, 0.7, 0.5)
  expected_hall_symbol = ' I 4'

  result = run_one_indexing(" ".join(pickle_paths),  " ".join(sweep_paths),
                            extra_args, expected_unit_cell,
                            expected_rmsds, expected_hall_symbol)
  assert len(result.indexed_reflections) > 1250, len(result.indexed_reflections)

def exercise_17():
  # test for small molecule multi-sweep indexing, 3 sweeps with different values
  # of goniometer setting rotation (i.e. phi scans)
  data_dir = os.path.join(dials_regression, "dials-191")
  import glob
  pickle_paths = [
    glob.glob(os.path.join(data_dir, "*SWEEP%i*_strong.pickle" %(i+1)))[0]
    for i in range(3)]
  sweep_paths = [
    glob.glob(os.path.join(data_dir, "*SWEEP%i*_datablock.json" %(i+1)))[0]
    for i in range(3)]
  extra_args = ["filter_ice=False"]
  expected_unit_cell = uctbx.unit_cell(
    (9.385, 15.219, 17.019, 90.076, 89.980, 100.816))
  expected_rmsds = (0.30, 0.32, 0.005 )
  expected_hall_symbol = ' P 1'

  result = run_one_indexing(" ".join(pickle_paths),  " ".join(sweep_paths),
                            extra_args, expected_unit_cell,
                            expected_rmsds, expected_hall_symbol)
  assert len(result.indexed_reflections) > 12000, len(result.indexed_reflections)
  # expect at least indexed 2000 reflections per experiment
  for i in range(3):
    assert (result.indexed_reflections['id'] == i).count(True) > 2000

def run(args):
  if not libtbx.env.has_module("dials_regression"):
    print "Skipping exercise_index_3D_FFT_simple: dials_regression not present"
    return

  exercises = (exercise_1, exercise_2, exercise_3, exercise_4, exercise_5,
               exercise_6, exercise_7, exercise_8, exercise_9, exercise_10,
               exercise_11, exercise_12, exercise_13, exercise_14, exercise_15,
               exercise_16, exercise_17)
  if len(args):
    args = [int(arg) for arg in args]
    for arg in args: assert arg > 0
    exercises = [exercises[arg-1] for arg in args]

  from libtbx import easy_mp

  nproc = easy_mp.get_processes(libtbx.Auto)
  nproc = min(nproc, len(exercises))

  def run_parallel(args):
    assert len(args) == 1
    exercise = args[0]
    exercise()

  easy_mp.parallel_map(
    func=run_parallel,
    iterable=[(e,) for e in exercises],
    processes=nproc)

if __name__ == '__main__':
  import sys
  from libtbx.utils import show_times_at_exit
  show_times_at_exit()
  run(sys.argv[1:])
