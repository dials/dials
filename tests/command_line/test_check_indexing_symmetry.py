from __future__ import annotations

from pathlib import Path

import procrunner


def test(dials_regression, tmp_path):
    experiments_path = Path(dials_regression) / "misc_test_data" / "i04-indexed.json"
    pickle_path = Path(dials_regression) / "misc_test_data" / "i04-indexed.pickle"

    assert experiments_path.is_file()
    assert pickle_path.is_file()

    result = procrunner.run(
        [
            "dials.check_indexing_symmetry",
            experiments_path,
            pickle_path,
            "d_min=4",
            "d_max=10",
            "symop_threshold=0.6",
        ],
        working_directory=tmp_path,
    )
    assert not result.returncode and not result.stderr
