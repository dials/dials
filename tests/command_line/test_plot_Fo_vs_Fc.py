from __future__ import annotations

import datetime

import procrunner
import pytest

# Conflicts with downstream dependencies causing warnings. Temporarily
# ignore this test failing. Remove after 2022-11-01
if datetime.date.today() < datetime.date(2022, 11, 1):
    timed_xfail = pytest.mark.xfail(
        reason="Downstream dependencies causing warnings on some platforms. Give some time for resolution."
    )
else:
    timed_xfail = lambda x: x


@timed_xfail
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
