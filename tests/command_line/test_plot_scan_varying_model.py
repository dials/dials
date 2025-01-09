from __future__ import annotations

from dials.command_line import plot_scan_varying_model


def test(dials_data, tmp_path, capsys):
    data_dir = dials_data("refinement_test_data", pathlib=True)
    plot_scan_varying_model.run(
        [
            str(data_dir / "glucose_isomerase_sv_refined.json"),
            f"output.directory={tmp_path}",
        ]
    )
    captured = capsys.readouterr()
    assert not captured.err
    output_dir = tmp_path / "scan-varying_model"
    assert (output_dir / "orientation.png").is_file()
    assert (output_dir / "unit_cell.png").is_file()
