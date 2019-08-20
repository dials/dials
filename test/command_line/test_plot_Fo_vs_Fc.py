from __future__ import absolute_import, division, print_function

import os
import procrunner


def test(dials_data, run_in_tmpdir):
    mtz_file = dials_data("dials_for_ed").join("refmac_final.mtz")
    result = procrunner.run(["dials.plot_Fo_vs_Fc", "hklin={0}".format(mtz_file)])

    assert os.path.isfile("Fo_vs_Fc.pdf")
    assert "|Fe| = 42.0" in result["stdout"].decode()
    assert not result.returncode and not result.stderr
