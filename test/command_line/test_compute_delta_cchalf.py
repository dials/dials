"""Tests for dials.compute_delta_cchalf."""
from __future__ import absolute_import, division, print_function

import os
import procrunner
import pytest
from dials.command_line.compute_delta_cchalf import phil_scope, CCHalfFromMTZ


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
    assert tmpdir.join("compute_delta_cchalf.html").check()


def test_compute_delta_cchalf_scaled_data_mtz(dials_data, tmpdir):
    """Test dials.compute_delta_cchalf on scaled data."""
    location = dials_data("l_cysteine_4_sweeps_scaled")
    refls = location.join("scaled_20_25.refl").strpath
    expts = location.join("scaled_20_25.expt").strpath

    # First export the data
    command = ["dials.export", refls, expts, "partiality_threshold=0.99"]
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
    assert tmpdir.join("compute_delta_cchalf.html").check()


def test_compute_delta_cchalf(dials_regression):
    """Test compute delta cchalf on an integrated mtz."""

    filename = os.path.join(
        dials_regression, "delta_cchalf_test_data", "test.XDS_ASCII.mtz"
    )
    params = phil_scope.extract()
    params.nbins = 1

    script = CCHalfFromMTZ(params, filename)

    assert script.statistics.mean_cchalf == pytest.approx(0.9458229988198384)
    assert list(script.statistics.cchalf_i) == pytest.approx(
        [0.7958683828042507, 0.9423789245583555]
    )
