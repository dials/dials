from __future__ import annotations

import procrunner
import pytest


@pytest.mark.parametrize("filename", ["image_15799_master.h5", "image_15799.nxs"])
def test_convert_to_cbf(dials_data, filename, tmp_path):
    result = procrunner.run(
        [
            "dials.import",
            dials_data("vmxi_thaumatin", pathlib=True) / filename,
            "image_range=1,10",
        ],
        working_directory=tmp_path,
    )
    result.check_returncode()
    assert not result.stderr
    assert tmp_path.joinpath("imported.expt").is_file()

    result = procrunner.run(
        ["dials.convert_to_cbf", "imported.expt"], working_directory=tmp_path
    )
    result.check_returncode()
    assert not result.stderr

    assert len(list(tmp_path.glob("as_cbf_*.cbf"))) == 10
