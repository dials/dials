from __future__ import absolute_import, division, print_function

import os
import pytest

from dials.util.options import OptionParser
from dials.util.options import flatten_datablocks

pytestmark = pytest.mark.skipif(
  not os.access('/dls/i04/data/2019/cm23004-1/20190109/Eiger', os.R_OK),
  reason='Test images not available')

@pytest.mark.xfail
def test_not_master_h5():
  data_h5 = '/dls/i04/data/2019/cm23004-1/20190109/Eiger/gw/Thaum/Thau_4/Thau_4_1_000001.h5'
  parser = OptionParser(read_datablocks=True, read_datablocks_from_images=True)
  params, options = parser.parse_args([data_h5])
  datablocks = flatten_datablocks(params.input.datablock)
  assert len(datablocks) == 0
