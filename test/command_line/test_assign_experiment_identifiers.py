"""
Test for dev.dials.assign_experiment_identifiers
"""
from __future__ import absolute_import, division, print_function

import os
import pytest
from libtbx import easy_run
from dials.array_family import flex
from dxtbx.serialize import load

class run_assign_identifiers(object):
  def __init__(self, pickle_path_list, sweep_path_list, extra_args):
    args = ["dev.dials.assign_experiment_identifiers"] + pickle_path_list + \
      sweep_path_list + extra_args
    command = " ".join(args)
    print(command)
    _ = easy_run.fully_buffered(command=command).raise_if_errors()
    assert os.path.exists("assigned_experiments.json")
    assert os.path.exists("assigned_reflections.pickle")

def test_assign_identifiers(dials_regression, run_in_tmpdir):
  """Test for dev.dials.assign_experiment_identifiers"""
  pickle_path_list = []
  sweep_path_list = []
  data_dir = os.path.join(dials_regression, "xia2-28")
  for i in [20, 25]:
    pickle_path_list.append(os.path.join(
      data_dir, str(i)+"_integrated.pickle"))
    sweep_path_list.append(os.path.join(
      data_dir, str(i)+"_integrated_experiments.json"))

  run_assign_identifiers(pickle_path_list, sweep_path_list, extra_args=[])

  r = flex.reflection_table.from_pickle("assigned_reflections.pickle")
  e = load.experiment_list("assigned_experiments.json", check_format=False)
  r.assert_experiment_identifiers_are_consistent(e)
  assert list(r.experiment_identifiers().values()) == ['0', '1']
  assert list(r.experiment_identifiers().keys()) == [0, 1]
  assert list(e.identifiers()) == ['0', '1']

  # now run again, with already assigned data
  pickle_path_list = ["assigned_reflections.pickle"]
  sweep_path_list = ["assigned_experiments.json"]
  run_assign_identifiers(pickle_path_list, sweep_path_list, extra_args=[])

  r = flex.reflection_table.from_pickle("assigned_reflections.pickle")
  e = load.experiment_list("assigned_experiments.json", check_format=False)
  r.assert_experiment_identifiers_are_consistent(e)
  assert list(r.experiment_identifiers().values()) == ['0', '1']
  assert list(r.experiment_identifiers().keys()) == [0, 1]
  assert list(e.identifiers()) == ['0', '1']

  # now run again, with adding more data
  pickle_path_list = ["assigned_reflections.pickle"]
  sweep_path_list = ["assigned_experiments.json"]
  for i in [30, 35]:
    pickle_path_list.append(os.path.join(data_dir,
      str(i)+"_integrated.pickle"))
    sweep_path_list.append(os.path.join(data_dir,
      str(i)+"_integrated_experiments.json"))

  run_assign_identifiers(pickle_path_list, sweep_path_list, extra_args=['identifiers="0 5 10 15"'])

  r = flex.reflection_table.from_pickle("assigned_reflections.pickle")
  e = load.experiment_list("assigned_experiments.json", check_format=False)
  r.assert_experiment_identifiers_are_consistent(e)
  assert list(r.experiment_identifiers().values()) == ['0', '5', '10', '15']
  assert list(r.experiment_identifiers().keys()) == [0, 1, 2, 3]
  assert list(e.identifiers()) == ['0', '5', '10', '15']
