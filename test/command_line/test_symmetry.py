from __future__ import absolute_import, division, print_function

import os

import procrunner


def test_symmetry(dials_regression, tmpdir):
    """Simple test to check that dials.symmetry completes"""

    result = procrunner.run(
        [
            "dials.symmetry",
            os.path.join(dials_regression, "xia2-28", "20_integrated_experiments.json"),
            os.path.join(dials_regression, "xia2-28", "20_integrated.pickle"),
            os.path.join(dials_regression, "xia2-28", "25_integrated_experiments.json"),
            os.path.join(dials_regression, "xia2-28", "25_integrated.pickle"),
        ],
        working_directory=tmpdir,
    )
    assert not result.returncode and not result.stderr
    assert tmpdir.join("symmetrized.refl").check()
    assert tmpdir.join("symmetrized.expt").check()
