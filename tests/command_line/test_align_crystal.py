from __future__ import annotations

import shutil
import subprocess


def test_align_crystal(dials_data, tmp_path):
    path = str(dials_data("experiment_test_data", pathlib=True))
    result = subprocess.run(
        [shutil.which("dials.align_crystal"), f"{path}/kappa_experiments.json"],
        cwd=tmp_path,
        capture_output=True,
    )
    assert not result.returncode and not result.stderr
    assert result.stdout.replace(b"\r\n", b"\n").endswith(
        b"""\
Angles between reciprocal cell axes and principal experimental axes:
+---------------------+--------+--------+--------+
| Experimental axis   |     a* |     b* |     c* |
|---------------------+--------+--------+--------|
| GON_PHI             | 88.712 | 88.317 |  1.76  |
| Beam                | 13.904 | 46.197 | 88.481 |
| GON_OMEGA           | 88.712 | 88.317 |  1.76  |
+---------------------+--------+--------+--------+

Angles between unit cell axes and principal experimental axes:
+---------------------+--------+--------+--------+
| Experimental axis   |      a |      b |      c |
|---------------------+--------+--------+--------|
| GON_PHI             | 89.484 | 88.801 |  1.76  |
| Beam                | 43.844 | 76.182 | 88.481 |
| GON_OMEGA           | 89.484 | 88.801 |  1.76  |
+---------------------+--------+--------+--------+

Independent solutions:
+----------------+------------------+-------------+-----------+
| Primary axis   | Secondary axis   |   GON_KAPPA |   GON_PHI |
|----------------+------------------+-------------+-----------|
| c* (6-fold)    | a* (2-fold)      |       4.324 |   -77.874 |
| c* (6-fold)    | a* (2-fold)      |      -4.324 |   106.075 |
+----------------+------------------+-------------+-----------+
"""
    )
