from __future__ import absolute_import, division, print_function

import six.moves.cPickle as pickle
import os
import procrunner

def test_scan_varying_refinement_of_a_multiple_panel_detector(dials_regression, tmpdir):
  from dials.array_family import flex

  tmpdir.chdir()

  result = procrunner.run_process([
      "dials.refine",
      os.path.join(dials_regression, "refinement_test_data", "i23_as_24_panel_barrel", 'experiments.json'),
      os.path.join(dials_regression, "refinement_test_data", "i23_as_24_panel_barrel", 'indexed.pickle'),
      "scan_varying=true",
      "history=history.pickle",
      "outlier.separate_blocks=False",
  ])
  assert result['exitcode'] == 0
  assert result['stderr'] == ''

  # there are plenty of things we could do with the refinement history, but
  # here just check that final RMSDs are low enough
  with open('history.pickle', 'rb') as f:
    history = pickle.load(f)
  final_rmsd = history['rmsd'][-1]
  assert final_rmsd[0] < 0.05
  assert final_rmsd[1] < 0.04
  assert final_rmsd[2] < 0.0002

  # also check that the used_in_refinement flag got set correctly
  rt = flex.reflection_table.from_pickle('refined.pickle')
  uir = rt.get_flags(rt.flags.used_in_refinement)
  assert uir.count(True) == history['num_reflections'][-1]
