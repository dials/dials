from __future__ import absolute_import, division, print_function

import os
import procrunner

def test(dials_regression, tmpdir):
  tmpdir.chdir()

  experiments_path = os.path.join(dials_regression, "misc_test_data", "i04-indexed.json")
  pickle_path = os.path.join(dials_regression, "misc_test_data", "i04-indexed.pickle")

  for path in (experiments_path, pickle_path):
    assert os.path.exists(path)

  result = procrunner.run_process([
      'dials.check_indexing_symmetry',
      experiments_path,
      pickle_path,
      'd_min=4',
      'd_max=10',
      'symop_threshold=0.6',
  ])
  assert result['exitcode'] == 0
  assert result['stderr'] == ''
