# test running data reduction programs on ssx data
from __future__ import annotations

import procrunner


def test_ssx_reduction(dials_data, tmp_path):
    """
    Check that dials.cosym, dials.scale, dials.export and dials.merge run
    successfully on ssx data.
    Also test a few smaller analysis programs.
    """
    ssx = dials_data("cunir_serial_processed", pathlib=True)
    refls = ssx / "integrated.refl"
    expts = ssx / "integrated.expt"

    result = procrunner.run(
        ["dials.cosym", expts, refls],
        working_directory=tmp_path,
    )
    assert not result.returncode and not result.stderr
    cosym_expts = tmp_path / "symmetrized.expt"
    cosym_refls = tmp_path / "symmetrized.refl"
    assert cosym_expts.is_file()
    assert cosym_refls.is_file()
    assert (tmp_path / "dials.cosym.html").is_file()

    result = procrunner.run(
        ["dials.scale", cosym_expts, cosym_refls],
        working_directory=tmp_path,
    )
    assert not result.returncode and not result.stderr
    scale_expts = tmp_path / "scaled.expt"
    scale_refls = tmp_path / "scaled.refl"
    assert scale_expts.is_file()
    assert scale_refls.is_file()
    assert (tmp_path / "dials.scale.html").is_file()

    result = procrunner.run(
        ["dials.export", scale_expts, scale_refls],
        working_directory=tmp_path,
    )
    assert not result.returncode and not result.stderr
    assert (tmp_path / "scaled.mtz").is_file()

    result = procrunner.run(
        ["dials.merge", scale_expts, scale_refls],
        working_directory=tmp_path,
    )
    assert not result.returncode and not result.stderr
    assert (tmp_path / "merged.mtz").is_file()

    result = procrunner.run(
        ["dials.damage_analysis", scale_expts, scale_refls],
        working_directory=tmp_path,
    )
    assert not result.returncode and not result.stderr
    assert (tmp_path / "dials.damage_analysis.html").is_file()

    result = procrunner.run(
        ["dials.compute_delta_cchalf", scale_expts, scale_refls],
        working_directory=tmp_path,
    )
    assert not result.returncode and not result.stderr
    assert (tmp_path / "compute_delta_cchalf.html").is_file()
