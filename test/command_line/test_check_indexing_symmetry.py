from __future__ import absolute_import, division, print_function

import os
import procrunner


def test(dials_regression, run_in_tmpdir):
    experiments_path = os.path.join(
        dials_regression, "misc_test_data", "i04-indexed.json"
    )
    pickle_path = os.path.join(dials_regression, "misc_test_data", "i04-indexed.pickle")

    for path in (experiments_path, pickle_path):
        assert os.path.exists(path)

    result = procrunner.run(
        [
            "dials.check_indexing_symmetry",
            experiments_path,
            pickle_path,
            "d_min=4",
            "d_max=10",
            "symop_threshold=0.6",
        ]
    )
    assert not result.returncode and not result.stderr
