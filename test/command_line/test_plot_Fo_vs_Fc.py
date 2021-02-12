import procrunner


def test(dials_data, tmp_path):
    mtz_file = dials_data("lysozyme_electron_diffraction").join("refmac_final.mtz")
    result = procrunner.run(
        ["dials.plot_Fo_vs_Fc", "hklin=" + mtz_file.strpath], working_directory=tmp_path
    )
    assert not result.returncode and not result.stderr
    assert tmp_path.joinpath("Fo_vs_Fc.pdf").is_file()
    assert "|Fe| = 42.0" in result.stdout.decode()
