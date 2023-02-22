from __future__ import annotations

import os
import shutil
import subprocess

from dials.command_line import plot_scan_varying_model


def test(dials_regression, tmp_path, capsys):
    plot_scan_varying_model.run(
        [
            os.path.join(
                dials_regression,
                "refinement_test_data",
                "multi_sweep_one_sample",
                "glucose_isomerase",
                "SWEEP1",
                "index",
                "sv_refined_experiments.json",
            ),
            f"output.directory={tmp_path}",
        ]
    )
    captured = capsys.readouterr()
    assert not captured.err
    output_dir = tmp_path / "scan-varying_model"
    assert output_dir.is_file("orientation.png")
    assert output_dir.is_file("unit_cell.png")
