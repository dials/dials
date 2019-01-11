from __future__ import absolute_import, division, print_function

import os
import pytest
import procrunner

@pytest.mark.parametrize("space_group", [None, "P 1"])
def test_cosym(regression_data, run_in_tmpdir, space_group):

  reg_path = regression_data("multi_crystal_proteinase_k").strpath

  command = ['dials.cosym', 'space_group='+str(space_group)]
  for i in [1, 2, 3, 4, 5, 7, 8, 10]:
    command.append(os.path.join(reg_path, "experiments_"+str(i)+".json"))
    command.append(os.path.join(reg_path, "reflections_"+str(i)+".pickle"))

  result = procrunner.run_process(command)
  assert result['exitcode'] == 0
  assert result['stderr'] == ''
  assert os.path.exists("reindexed_reflections.pickle")
  assert os.path.exists("reindexed_experiments.json")
