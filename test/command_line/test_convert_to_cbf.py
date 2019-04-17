from __future__ import absolute_import, division, print_function

import procrunner
import pytest


@pytest.mark.parametrize("filename", ["image_15799_master.h5", "image_15799.nxs"])
def test_convert_to_cbf(dials_data, filename, tmpdir):
    result = procrunner.run(
        ["dials.import", dials_data("vmxi_thaumatin") / filename, "image_range=1,10"],
        working_directory=tmpdir,
    )
    result.check_returncode()
    assert not result.stderr
    assert tmpdir.join("imported_experiments.json").check()

    result = procrunner.run(
        ["dials.convert_to_cbf", "imported_experiments.json"], working_directory=tmpdir
    )
    result.check_returncode()
    assert not result.stderr

    assert len(tmpdir.listdir("as_cbf_*.cbf")) == 10
