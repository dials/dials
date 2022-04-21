from __future__ import annotations

import os

import procrunner
import pytest

cpp_tests = [
    # Paths are under /build/
    "tests/algorithms/spatial_indexing/tst_collision_detection",
    "tests/algorithms/spot_prediction/tst_reeke_model",
]


@pytest.mark.parametrize(
    "executable", cpp_tests, ids=[p.replace("/", "-") for p in cpp_tests]
)
def test_cpp_program(executable):
    full_path = os.path.join(
        os.environ["LIBTBX_BUILD"], "dials", *(executable.split("/"))
    )
    print(full_path)

    result = procrunner.run([full_path])
    assert not result.returncode and not result.stderr
