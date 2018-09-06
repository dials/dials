from __future__ import absolute_import, division, print_function

import six.moves.cPickle as pickle
import os

import procrunner
import pytest

def test_joint_refinement(dials_regression, tmpdir):
  """A basic test of joint refinement of the CS-PAD detector at hierarchy level 2
     with 300 crystals."""
  from dials.array_family import flex

  bevington = pytest.importorskip("scitbx.examples.bevington")
  if not hasattr(bevington, 'non_linear_ls_eigen_wrapper'):
    pytest.skip("Skipping test as SparseLevMar engine not available")

  tmpdir.chdir()

  data_dir = os.path.join(dials_regression, "refinement_test_data", "xfel_metrology")

  # Do refinement and load the history
  result = procrunner.run_process([
      "dials.refine",
      os.path.join(data_dir, 'benchmark_level2d.json'),
      os.path.join(data_dir, 'benchmark_level2d.pickle'),
      os.path.join(data_dir, 'refine.phil'),
      "history=history.pickle",
  ])
  assert result['exitcode'] == 0
  assert result['stderr'] == ''

  # there are plenty of things we could do with the refinement history, but
  # here just check that final RMSDs are low enough
  with open('history.pickle', 'rb') as f:
    history = pickle.load(f)
  final_rmsd = history['rmsd'][-1]
  assert final_rmsd[0] < 0.0354
  assert final_rmsd[1] < 0.0406
  assert final_rmsd[2] < 0.0018

  # also check that the used_in_refinement flag got set correctly
  rt = flex.reflection_table.from_pickle('refined.pickle')
  uir = rt.get_flags(rt.flags.used_in_refinement)
  assert uir.count(True) == history['num_reflections'][-1]
