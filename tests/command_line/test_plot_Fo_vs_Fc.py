from __future__ import annotations

from dials.command_line import plot_Fo_vs_Fc


def test(dials_data, tmp_path, capsys):
    mtz_file = dials_data("lysozyme_electron_diffraction") / "refmac_final.mtz"
    plot_file = tmp_path / "Fo_vs_Fc.pdf"
    plot_Fo_vs_Fc.run(
        [f"hklin={mtz_file}", f"plot_filename={plot_file}"],
    )
    captured = capsys.readouterr()
    assert not captured.err
    assert plot_file.is_file()
    assert "|Fe| = 42.0" in captured.out
