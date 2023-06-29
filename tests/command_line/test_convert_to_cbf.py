from __future__ import annotations

import os
import shutil
import subprocess

import pytest


@pytest.mark.parametrize("filename", ["image_15799_master.h5", "image_15799.nxs"])
def test_convert_to_cbf(dials_data, filename, tmp_path):
    result = subprocess.run(
        [
            shutil.which("dials.import"),
            dials_data("vmxi_thaumatin", pathlib=True) / filename,
            "image_range=1,10",
        ],
        cwd=tmp_path,
        capture_output=True,
    )
    result.check_returncode()
    assert not result.stderr
    assert tmp_path.joinpath("imported.expt").is_file()

    result = subprocess.run(
        [shutil.which("dials.convert_to_cbf"), "imported.expt"],
        cwd=tmp_path,
        capture_output=True,
        env={
            **os.environ,
            "PYTHONWARNINGS": ",".join([
                "ignore:`product` is deprecated as of NumPy 1.25.0:DeprecationWarning",
                "ignore:pkg_resources is deprecated as an API.:DeprecationWarning",
            ])
        },
    )
    result.check_returncode()
    assert not result.stderr

    assert len(list(tmp_path.glob("as_cbf_*.cbf"))) == 10
