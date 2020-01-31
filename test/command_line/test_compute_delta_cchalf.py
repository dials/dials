"""Tests for dials.compute_delta_cchalf."""
from __future__ import absolute_import, division, print_function

import procrunner


def check_cchalf_result(fileobj):
    """Inspect the result"""
    lines = fileobj.readlines()
    assert lines[0] == "1 -0.004673\n"
    assert lines[1] == "0 0.001234\n"


def test_compute_delta_cchalf_scaled_data(dials_data, tmpdir):
    """Test dials.compute_delta_cchalf on scaled data."""
    location = dials_data("l_cysteine_4_sweeps_scaled")
    refls = location.join("scaled_20_25.refl").strpath
    expts = location.join("scaled_20_25.expt").strpath

    # set cutoff to 0.0 to force one to be 'rejected'
    command = [
        "dials.compute_delta_cchalf",
        refls,
        expts,
        "stdcutoff=0.0",
        "output.reflections=filtered.refl",
        "output.experiments=filtered.expt",
    ]
    result = procrunner.run(command, working_directory=tmpdir)
    assert not result.returncode and not result.stderr
    assert tmpdir.join("filtered.expt").check()
    assert tmpdir.join("filtered.refl").check()
    assert tmpdir.join("delta_cchalf.dat").check()
    assert tmpdir.join("compute_delta_cchalf.html").check()
    with open(tmpdir.join("delta_cchalf.dat").strpath, "r") as f:
        check_cchalf_result(f)


def test_compute_delta_cchalf_scaled_data_mtz(dials_data, tmpdir):
    """Test dials.compute_delta_cchalf on scaled data."""
    location = dials_data("l_cysteine_4_sweeps_scaled")
    refls = location.join("scaled_20_25.refl").strpath
    expts = location.join("scaled_20_25.expt").strpath

    # First export the data
    command = ["dials.export", refls, expts]
    result = procrunner.run(command, working_directory=tmpdir)
    assert not result.returncode and not result.stderr
    assert tmpdir.join("scaled.mtz").check()

    # set cutoff to 0.0 to force one to be 'rejected'
    command = [
        "dials.compute_delta_cchalf",
        "mtzfile=%s" % tmpdir.join("scaled.mtz").strpath,
        "stdcutoff=0.0",
    ]
    result = procrunner.run(command, working_directory=tmpdir)
    assert not result.returncode and not result.stderr
    assert tmpdir.join("delta_cchalf.dat").check()
    assert tmpdir.join("compute_delta_cchalf.html").check()
    with open(tmpdir.join("delta_cchalf.dat").strpath, "r") as f:
        check_cchalf_result(f)
