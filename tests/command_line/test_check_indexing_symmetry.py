from __future__ import annotations

import shutil
import subprocess
from pathlib import Path


def test(dials_regression, tmp_path):
    experiments_path = Path(dials_regression) / "misc_test_data" / "i04-indexed.json"
    pickle_path = Path(dials_regression) / "misc_test_data" / "i04-indexed.pickle"

    assert experiments_path.is_file()
    assert pickle_path.is_file()

    result = subprocess.run(
        [
            shutil.which("dials.check_indexing_symmetry"),
            experiments_path,
            pickle_path,
            "d_min=4",
            "d_max=10",
            "symop_threshold=0.6",
        ],
        cwd=tmp_path,
        capture_output=True,
    )
    assert not result.returncode and not result.stderr
