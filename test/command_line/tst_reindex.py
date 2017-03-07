from __future__ import absolute_import, division
import os
import libtbx.load_env
from libtbx import easy_run
from libtbx import easy_pickle
from libtbx.test_utils import approx_equal
from libtbx.test_utils import open_tmp_directory
from dxtbx.serialize import load

have_dials_regression = libtbx.env.has_module("dials_regression")
if have_dials_regression:
  dials_regression = libtbx.env.find_in_repositories(
    relative_path="dials_regression",
    test=os.path.isdir)


def run():
  if not have_dials_regression:
    print "Skipping tst_reindex: dials_regression not available."
    return

  data_dir = os.path.join(dials_regression, "indexing_test_data", "i04_weak_data")
  pickle_path = os.path.join(data_dir, "indexed.pickle")
  experiments_path = os.path.join(data_dir, "experiments.json")
  commands = ["dials.reindex",
              pickle_path,
              experiments_path,
              "change_of_basis_op=2a,b,c",
              "space_group=P1"]
  command = " ".join(commands)
  print command
  cwd = os.path.abspath(os.curdir)
  tmp_dir = open_tmp_directory()
  os.chdir(tmp_dir)
  result = easy_run.fully_buffered(command=command).raise_if_errors()
  old_reflections = easy_pickle.load(pickle_path)
  assert os.path.exists('reindexed_reflections.pickle')
  new_reflections = easy_pickle.load('reindexed_reflections.pickle')
  old_experiments = load.experiment_list(experiments_path, check_format=False)
  assert os.path.exists('reindexed_experiments.json')
  new_experiments = load.experiment_list(
    'reindexed_experiments.json', check_format=False)
  h1,k1,l1 = old_reflections['miller_index'].as_vec3_double().parts()
  h2,k2,l2 = new_reflections['miller_index'].as_vec3_double().parts()
  assert approx_equal(2*h1, h2)
  assert approx_equal(k1, k2)
  assert approx_equal(l1, l2)
  old_uc_params = old_experiments[0].crystal.get_unit_cell().parameters()
  new_uc_params = new_experiments[0].crystal.get_unit_cell().parameters()
  assert approx_equal(2*old_uc_params[0], new_uc_params[0])
  assert approx_equal(old_uc_params[1:], new_uc_params[1:])
  assert old_experiments[0].crystal.get_space_group().type().hall_symbol() == ' P 1'
  assert new_experiments[0].crystal.get_space_group().type().hall_symbol() == ' P 1'

  # set space group P4
  from cctbx import sgtbx
  cb_op = sgtbx.change_of_basis_op('a,b,c')
  commands = ["dials.reindex",
              experiments_path,
              "space_group=P4",
              "change_of_basis_op=%s" %str(cb_op),
              "output.experiments=P4.json"]
  command = " ".join(commands)
  print command
  result = easy_run.fully_buffered(command=command).raise_if_errors()
  # apply one of the symops from the space group
  cb_op = sgtbx.change_of_basis_op('-x,-y,z')
  commands = ["dials.reindex",
              "P4.json",
              "change_of_basis_op=%s" %str(cb_op),
              "output.experiments=P4_reindexed.json"]
  command = " ".join(commands)
  print command
  result = easy_run.fully_buffered(command=command).raise_if_errors()
  new_experiments1 = load.experiment_list(
    'P4_reindexed.json', check_format=False)
  assert approx_equal(
    new_experiments1[0].crystal.get_A(),
    old_experiments[0].crystal.change_basis(cb_op).get_A())
  #
  cb_op = sgtbx.change_of_basis_op('-x,-y,z')
  commands = ["dials.reindex",
              "P4.json",
              "change_of_basis_op=auto",
              "reference=P4_reindexed.json",
              "output.experiments=P4_reindexed2.json"]
  command = " ".join(commands)
  print command
  result = easy_run.fully_buffered(command=command).raise_if_errors()
  new_experiments2 = load.experiment_list(
    'P4_reindexed2.json', check_format=False)
  assert approx_equal(
    new_experiments1[0].crystal.get_A(),
    new_experiments2[0].crystal.get_A())



if __name__ == '__main__':
  run()
  print "OK"
