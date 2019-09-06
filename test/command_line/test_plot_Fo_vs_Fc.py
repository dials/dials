from __future__ import absolute_import, division, print_function

import procrunner


def test(dials_data, tmpdir):
    mtz_file = dials_data("lysozyme_electron_diffraction").join("refmac_final.mtz")
    result = procrunner.run(
        ["dials.plot_Fo_vs_Fc", "hklin={0}".format(mtz_file.strpath)],
        working_directory=tmpdir,
    )
    assert not result.returncode and not result.stderr
    assert tmpdir.join("Fo_vs_Fc.pdf").check()
    assert "|Fe| = 42.0" in result["stdout"].decode()
