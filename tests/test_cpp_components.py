from __future__ import annotations

import os
import subprocess

import pytest

cpp_tests = [
    # Paths are under /build/
    "tests/algorithms/spatial_indexing/tst_collision_detection",
    "tests/algorithms/spot_prediction/tst_reeke_model",
    "tests/util/tst_thread_pool",
]


@pytest.mark.parametrize(
    "executable", cpp_tests, ids=[p.replace("/", "-") for p in cpp_tests]
)
def test_cpp_program(executable):
    if "LIBTBX_BUILD" not in os.environ:
        pytest.skip("LIBTBX_ENV is unset; don't know how to find test executable")
    full_path = os.path.join(
        os.environ["LIBTBX_BUILD"], "dials", *(executable.split("/"))
    )
    print(full_path)

    result = subprocess.run([full_path], capture_output=True)
    assert not result.returncode and not result.stderr
