from __future__ import absolute_import, division, print_function

import procrunner
import pytest


def test_rs_mapper(dials_data, tmpdir):
    result = procrunner.run(
        [
            "dials.rs_mapper",
            dials_data("centroid_test_data").join("datablock.json").strpath,
            'map_file="junk.ccp4"',
        ],
        working_directory=tmpdir.strpath,
    )
    assert not result.returncode and not result.stderr
    assert tmpdir.join("junk.ccp4").check()

    # load results
    from iotbx import ccp4_map
    from scitbx.array_family import flex

    m = ccp4_map.map_reader(file_name=tmpdir.join("junk.ccp4").strpath)
    assert len(m.data) == 7189057
    assert m.header_min == -1.0
    assert flex.min(m.data) == -1.0

    assert m.header_max == 2052.75
    assert flex.max(m.data) == 2052.75

    assert m.header_mean == pytest.approx(0.018606403842568398, abs=1e-6)
    assert flex.mean(m.data) == pytest.approx(0.018606403842568398, abs=1e-6)
