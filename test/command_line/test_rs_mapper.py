from __future__ import absolute_import, division, print_function

import os
import pytest
from libtbx import easy_run

def test_rs_mapper(dials_regression, tmpdir):
  tmpdir.chdir()

  data_dir = os.path.join(dials_regression, "centroid_test_data")
  datablock_path = os.path.join(data_dir, "datablock.json")

  cmd = 'dials.rs_mapper ' + datablock_path + ' map_file="junk.ccp4"'

  result = easy_run.fully_buffered(command=cmd).raise_if_errors()
  # load results
  from iotbx import ccp4_map
  from scitbx.array_family import flex
  m = ccp4_map.map_reader(file_name="junk.ccp4")
  assert len(m.data) == 7189057
  assert m.header_min == -1.0
  assert flex.min(m.data) == -1.0

  assert m.header_max == 2052.75
  assert flex.max(m.data) == 2052.75

  assert m.header_mean == pytest.approx(0.018606403842568398)
  assert flex.mean(m.data) == pytest.approx(0.018606403842568398)
