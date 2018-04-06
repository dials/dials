from __future__ import absolute_import, division, print_function

import os
import procrunner

def test(dials_regression, tmpdir):
  tmpdir.chdir()

  input_filename = os.path.join(dials_regression, "centroid_test_data", "datablock.json")

  result = procrunner.run_process([
      'dials.estimate_gain',
      'input.datablock=' + input_filename,
  ])
  assert result['exitcode'] == 0
  assert result['stderr'] == ''
  assert 'Estimated gain: 1.0' in result['stdout']
