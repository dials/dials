from __future__ import division
import os
import libtbx.load_env
from libtbx import easy_run
from libtbx.test_utils import approx_equal
from scitbx import matrix
from cctbx.array_family import flex # import dependency
from cctbx import uctbx
from dials.model.serialize import load


def exercise_index_3D_FFT_simple():

  if not libtbx.env.has_module("dials_regression"):
    print "Skipping exercise_index_3D_FFT_simple: dials_regression not present"
    return

  data_dir = libtbx.env.find_in_repositories(
    relative_path="dials_regression/centroid_test_data",
    test=os.path.isdir)
  template = os.path.join(data_dir, "centroid_00*.cbf")
  args = ["dials.spotfinder", "min_spot_size=2",
          template, "-o reflections.pickle"]
  result = easy_run.fully_buffered(command=" ".join(args)).raise_if_errors()
  assert os.path.exists("reflections.pickle")
  reflections = load.reflections("reflections.pickle")

  py_script_path = libtbx.env.find_in_repositories(
    relative_path="dials/scratch/rjg/index_3D_FFT_simple.py",
    test=os.path.isfile)
  args = ["dials.python", py_script_path, "reflections.pickle", template,
          "fix_detector=True", "fix_beam=True", "n_macro_cycles=10"]
  result = easy_run.fully_buffered(command=" ".join(args)).raise_if_errors()
  assert os.path.exists("indexed.pickle")
  reflections = load.reflections("indexed.pickle")
  mi = reflections.miller_index()
  assert (mi != (0,0,0)).count(True) >= 720
  assert (mi != (0,0,0)).count(False) <= 11
  assert os.path.exists("crystal.json")
  crystal_model = load.crystal("crystal.json")
  expected_unit_cell = uctbx.unit_cell(
    (39.5, 42.1, 42.1, 90, 90, 90))
  assert crystal_model.get_unit_cell().is_similar_to(
    expected_unit_cell,
    relative_length_tolerance=0.005,
    absolute_angle_tolerance=0.5)
  sg = crystal_model.get_space_group()
  assert sg.type().hall_symbol() == ' P 1'
  assert os.path.exists("sweep.json")
  sweep = load.sweep("sweep.json")
  reflections = reflections.select(mi != (0,0,0))
  rmsds = get_rmsds_obs_pred(reflections, sweep, crystal_model)
  expected_rmsds = (0.187, 0.174, 0.00353)
  for actual, expected in zip(rmsds, expected_rmsds):
    assert actual <= expected

  # now enforce symmetry
  args = ["dials.python", py_script_path, "reflections.pickle", template,
          "fix_detector=True", "fix_beam=True", "n_macro_cycles=10",
          "space_group=P4"]
  result = easy_run.fully_buffered(command=" ".join(args)).raise_if_errors()
  assert os.path.exists("indexed.pickle")
  reflections = load.reflections("indexed.pickle")
  mi = reflections.miller_index()
  assert (mi != (0,0,0)).count(True) >= 719
  assert (mi != (0,0,0)).count(False) <= 12
  assert os.path.exists("crystal.json")
  crystal_model = load.crystal("crystal.json")
  expected_unit_cell = uctbx.unit_cell(
    (42.1, 42.1, 39.5, 90, 90, 90))
  assert crystal_model.get_unit_cell().is_similar_to(
    expected_unit_cell,
    relative_length_tolerance=0.005,
    absolute_angle_tolerance=0.1)
  a, b, c = crystal_model.get_real_space_vectors()
  assert approx_equal(a.length(), b.length())
  assert c.length() < b.length()
  assert approx_equal(a.angle(b, deg=True), 90)
  assert approx_equal(b.angle(c, deg=True), 90)
  assert approx_equal(c.angle(a, deg=True), 90)
  sg = crystal_model.get_space_group()
  assert sg.type().hall_symbol() == ' P 4'
  assert os.path.exists("sweep.json")
  sweep = load.sweep("sweep.json")
  reflections = reflections.select(mi != (0,0,0))
  rmsds = get_rmsds_obs_pred(reflections, sweep, crystal_model)
  expected_rmsds = (0.210, 0.255, 0.00370)
  for actual, expected in zip(rmsds, expected_rmsds):
    assert actual <= expected

def get_rmsds_obs_pred(observations, sweep, crystal_model):
  from dials.algorithms.spot_prediction import ray_intersection
  reflections = ray_intersection(sweep.get_detector(), observations)
  from dials.scratch.rjg.index_3D_FFT_simple import master_params
  from dials.algorithms.refinement import RefinerFactory
  refine = RefinerFactory.from_parameters(master_params, verbosity=0)
  refine.prepare(sweep, crystal_model, reflections)
  return refine.rmsds()

def run():
  exercise_index_3D_FFT_simple()

if __name__ == '__main__':
  from libtbx.utils import show_times_at_exit
  show_times_at_exit()
  run()
