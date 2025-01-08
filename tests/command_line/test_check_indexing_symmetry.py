from __future__ import annotations

import shutil
import subprocess


def test(dials_data, tmp_path):
    data_dir = dials_data("insulin_processed", pathlib=True)
    experiments_path = data_dir / "indexed.expt"
    pickle_path = data_dir / "indexed.refl"

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
