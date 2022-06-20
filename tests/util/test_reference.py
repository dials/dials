from __future__ import annotations

import os

import pytest

from dials.util.reference import (
    intensities_from_reference_data_file,
    intensities_from_reference_file,
    intensities_from_reference_model_file,
)


def test_intensities_from_reference_model(dials_data):
    "Test importing from a pdb model"

    pdb_file = os.fspath(dials_data("cunir_serial", pathlib=True) / "2bw4.pdb")
    intensities = intensities_from_reference_file(pdb_file)
    assert intensities.data()
    assert not intensities.anomalous_flag()

    i2 = intensities_from_reference_model_file(pdb_file, wavelength=1.0)
    assert i2.data()
    assert i2.anomalous_flag()

    bad_input = "bad.pdbb"
    with pytest.raises(ValueError):
        _ = intensities_from_reference_file(bad_input)
    with pytest.raises(ValueError):
        _ = intensities_from_reference_model_file(bad_input)


def test_intensities_from_reference_data_file(dials_data):
    "Test importing from an mtz datafile"

    mtz_file = os.fspath(
        dials_data("x4wide_processed", pathlib=True) / "AUTOMATIC_DEFAULT_scaled.mtz"
    )
    intensities = intensities_from_reference_file(mtz_file)
    assert intensities.data()
    assert intensities.anomalous_flag()

    i2 = intensities_from_reference_data_file(mtz_file)
    assert i2.data() == intensities.data()

    bad_input = "bad.mtzz"
    with pytest.raises(ValueError):
        _ = intensities_from_reference_data_file(bad_input)
