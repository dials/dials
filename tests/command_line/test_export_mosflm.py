from __future__ import annotations

import shutil
import subprocess


def test_export_mosflm(dials_data, tmp_path):
    data_dir = dials_data("i04_weak_data")
    experiments = data_dir / "experiments.json"
    result = subprocess.run(
        [shutil.which("dials.export"), "format=mosflm", experiments],
        cwd=tmp_path,
        capture_output=True,
    )
    assert not result.returncode and not result.stderr

    assert (tmp_path / "mosflm" / "index.mat").is_file()
    lines = (tmp_path / "mosflm" / "index.mat").read_text()
    assert (
        lines
        == """
  0.01492625  0.00495607 -0.00238037
  0.00660784 -0.01505459  0.00149210
 -0.00437699 -0.00583805 -0.00587092
       0.000       0.000       0.000
  0.88323967  0.29347317 -0.36573377
  0.39109500 -0.89133994  0.22925490
 -0.25871295 -0.34552367 -0.90204268
     57.7658     57.7991    149.9967     90.0076     90.0138     90.0101
       0.000       0.000       0.000
""".strip("\n")
    )
    assert (tmp_path / "mosflm" / "mosflm.in").is_file()
    lines = (tmp_path / "mosflm" / "mosflm.in").read_text()
    assert (
        lines
        == """
TEMPLATE th_8_2_####.cbf
SYMMETRY 1
BEAM 205.383 211.023
DISTANCE 265.1161
MATRIX index.mat
""".strip("\n")
    )
