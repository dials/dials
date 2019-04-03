from __future__ import absolute_import, division, print_function

import glob
import os
import pytest

import procrunner

pytestmark = pytest.mark.skipif(
    not os.access("/dls/i04/data/2019/cm23004-1/20190109/Eiger", os.R_OK),
    reason="Test images not available",
)


@pytest.mark.parametrize(
    "master_h5",
    [
        "/dls/i04/data/2019/cm23004-1/20190109/Eiger/gw/Thaum/Thau_4/Thau_4_1_master.h5",
        "/dls/i04/data/2019/cm23004-1/20190109/Eiger/gw/Thaum/Thau_4/Thau_4_1.nxs",
    ],
)
def test_convert_to_cbf(master_h5):
    result = procrunner.run(["dials.convert_to_cbf", master_h5])
    assert result["exitcode"] == 0
    assert result["stderr"] == ""
    g = glob.glob("as_cbf_*.cbf")
    assert len(g) == 900  # need a smaller test set!
