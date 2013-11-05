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



class run_one_indexing(object):
  def __init__(self,
               pickle_path,
               sweep_path,
               extra_args,
               expected_unit_cell,
               expected_rmsds,
               expected_hall_symbol):

    py_script_path = libtbx.env.find_in_repositories(
      relative_path="dials/scratch/rjg/index_3D_FFT_simple.py",
      test=os.path.isfile)

    args = ["dials.python", py_script_path, pickle_path, sweep_path] + extra_args

    cwd = os.path.abspath(os.curdir)
    tmp_dir = open_tmp_directory(suffix="test_3DFFT_indexing")
    os.chdir(tmp_dir)
    result = easy_run.fully_buffered(command=" ".join(args)).raise_if_errors()
    os.chdir(cwd)
    assert os.path.exists(os.path.join(tmp_dir, "indexed.pickle"))
    assert os.path.exists(os.path.join(tmp_dir, "crystal.json"))
    assert os.path.exists(os.path.join(tmp_dir, "sweep.json"))
    self.reflections = load.reflections(os.path.join(tmp_dir, "indexed.pickle"))
    self.crystal_model = load.crystal(os.path.join(tmp_dir, "crystal.json"))
    self.sweep = load.sweep(os.path.join(tmp_dir, "sweep.json"))
    assert self.crystal_model.get_unit_cell().is_similar_to(
      expected_unit_cell,
      relative_length_tolerance=0.005,
      absolute_angle_tolerance=0.5)
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
    refine = RefinerFactory.from_parameters(master_params, verbosity=0)
    refine.prepare(sweep, crystal_model, reflections)
    return refine.rmsds()

def exercise_index_3D_FFT_simple():

  if not libtbx.env.has_module("dials_regression"):
    print "Skipping exercise_index_3D_FFT_simple: dials_regression not present"
    return

  dials_regression = libtbx.env.find_in_repositories(
    relative_path="dials_regression",
    test=os.path.isdir)

  # thaumatin
  data_dir = os.path.join(dials_regression, "indexing_test_data", "i04-weak-data")
  pickle_path = os.path.join(data_dir, "full.pickle")
  sweep_path = os.path.join(data_dir, "sweep_orig.json")
  extra_args = ["multiple_lattice_search=False", # use older non-clustering version
                "reflections_per_degree=5"]
  expected_unit_cell = uctbx.unit_cell(
    (58, 58, 150, 90, 90, 90))
  expected_rmsds = (0.06, 0.05, 0.0004)
  expected_hall_symbol = ' P 1'

  result = run_one_indexing(pickle_path, sweep_path, extra_args, expected_unit_cell,
                     expected_rmsds, expected_hall_symbol)

  extra_args = ["multiple_lattice_search=True",
                "reflections_per_degree=5",
                "d_min=4"]

  result = run_one_indexing(pickle_path, sweep_path, extra_args, expected_unit_cell,
                            expected_rmsds, expected_hall_symbol)

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


def run():
  exercise_index_3D_FFT_simple()

if __name__ == '__main__':
  from libtbx.utils import show_times_at_exit
  show_times_at_exit()
  run()
