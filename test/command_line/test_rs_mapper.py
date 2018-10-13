from __future__ import absolute_import, division, print_function

import os

import procrunner
import pytest

def test_rs_mapper(dials_regression, run_in_tmpdir):
  result = procrunner.run_process([
      'dials.rs_mapper',
      os.path.join(dials_regression, "centroid_test_data", "datablock.json"),
      'map_file="junk.ccp4"',
  ])
  assert result['exitcode'] == 0
  assert result['stderr'] == ''
  assert os.path.exists('junk.ccp4')

  # load results
  from iotbx import ccp4_map
  from scitbx.array_family import flex
  m = ccp4_map.map_reader(file_name="junk.ccp4")
  assert len(m.data) == 7189057
  assert m.header_min == -1.0
  assert flex.min(m.data) == -1.0

  assert m.header_max == 2052.75
  assert flex.max(m.data) == 2052.75

  assert m.header_mean == pytest.approx(0.018606403842568398, abs=1e-6)
  assert flex.mean(m.data) == pytest.approx(0.018606403842568398, abs=1e-6)
