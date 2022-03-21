"""
Tests for dials.command_line.anvil_correction.
"""


from __future__ import annotations

import copy
from pathlib import Path

import pytest

from dxtbx.model import ExperimentList

from dials.array_family import flex
from dials.command_line.anvil_correction import (
    correct_intensities_for_dac_attenuation,
    run,
)


def test_correct_correction(dials_data):
    """Test that the anvil absorption correction is producing expected values."""
    data_dir = dials_data("centroid_test_data", pathlib=True)

    # We'll need an integrated reflection table and an experiment list.
    reflections_file = data_dir / "integrated.pickle"
    experiments_file = data_dir / "experiments.json"

    # We need only test with the first ten reflections.
    reflections = flex.reflection_table.from_file(reflections_file)
    reflections = reflections.select(flex.size_t_range(10))

    experiment = ExperimentList.from_file(experiments_file)[0]

    # Test the correction that would be applied to a DAC with 1.5mm-thick anvils,
    # aligned along the z-axis at goniometer zero-datum.
    old_reflections = copy.deepcopy(reflections)
    correct_intensities_for_dac_attenuation(experiment, reflections, (0, 0, 1), 1.5)

    cases = {
        "intensity.sum.value": reflections.flags.integrated_sum,
        "intensity.sum.variance": reflections.flags.integrated_sum,
        "intensity.prf.value": reflections.flags.integrated_prf,
        "intensity.prf.variance": reflections.flags.integrated_prf,
    }
    corrections = flex.double(
        [
            0,
            6.653068275094517,
            6.522657529202368,
            6.3865190053761,
            6.587270967838122,
            6.43403642876391,
            6.39216742203502,
            0,
            6.152148372872684,
            6.0474840161407375,
        ]
    )
    for case, flag in cases.items():
        flagged = reflections.get_flags(flag)

        target_correction = corrections.select(flagged)
        if "variance" in case:
            target_correction = flex.pow2(target_correction)

        intensity_correction = (reflections[case] / old_reflections[case]).select(
            flagged
        )

        # Check that the un-integrated reflections are unchanged.
        assert pytest.approx(reflections[case].select(~flagged)) == old_reflections[
            case
        ].select(~flagged), (
            "Un-integrated reflections have been erroneously " "'corrected'."
        )

        # Check that the applied corrections are correct.
        assert pytest.approx(intensity_correction, rel=1e-5) == list(
            target_correction
        ), f"The applied intensity correction to {case} doesn't seem to be correct."


def test_help_message(dials_data, capsys):
    """Test that we get a help message when improper input is provided."""
    data_dir = dials_data("centroid_test_data", pathlib=True)

    # We'll need an integrated reflection table and an experiment list.
    reflections_file = str(data_dir / "integrated.pickle")
    experiments_file = str(data_dir / "experiments.json")

    for arguments in (
        [],
        [reflections_file],
        [experiments_file],
        [experiments_file, reflections_file, "anvil.normal=0,0,0"],
    ):
        with pytest.raises(SystemExit):
            run(arguments)
            assert (
                "Usage: dials.anvil_correction [options] integrated.expt "
                "integrated.refl" in capsys.readouterr().err
            )


def test_command_line(dials_data, run_in_tmp_path):
    """Test that the program runs correctly."""
    data_dir = dials_data("centroid_test_data", pathlib=True)

    # We'll need an integrated reflection table and an experiment list.
    reflections_file = str(data_dir / "integrated.pickle")
    experiments_file = str(data_dir / "experiments.json")

    run([experiments_file, reflections_file])

    assert Path("corrected.refl").is_file()

    with Path("dials.anvil_correction.log").open() as f:
        logfile = f.read()

    assert "Correcting integrated reflection intensities" in logfile
    assert "Writing the reflection table" in logfile
