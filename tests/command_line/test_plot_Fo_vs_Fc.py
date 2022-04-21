from __future__ import annotations

import procrunner


def test(dials_data, tmp_path):
    mtz_file = (
        dials_data("lysozyme_electron_diffraction", pathlib=True) / "refmac_final.mtz"
    )
    result = procrunner.run(
        ["dials.plot_Fo_vs_Fc", f"hklin={mtz_file}"], working_directory=tmp_path
    )
    assert not result.returncode and not result.stderr
    assert tmp_path.joinpath("Fo_vs_Fc.pdf").is_file()
    assert "|Fe| = 42.0" in result.stdout.decode()
