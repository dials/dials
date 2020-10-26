"""
Tests for dials.command_line.anvil_correction.
"""

from __future__ import absolute_import, division, print_function

import copy

import pytest

from dxtbx.model import ExperimentList

from dials.array_family import flex
from dials.command_line.anvil_correction import (
    correct_intensities_for_dac_attenuation,
    run,
)


def test_correct_correction(dials_data):
    """Test that the anvil absorption correction is producing expected values."""
    data_dir = dials_data("centroid_test_data")

    # We'll need an integrated reflection table and an experiment list.
    reflections_file = data_dir.join("integrated.pickle")
    experiments_file = data_dir.join("experiments.json")

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
        ), ("The applied intensity correction to %s doesn't seem to be correct." % case)


def test_help_message(dials_data, capsys):
    """Test that we get a help message when improper input is provided."""
    data_dir = dials_data("centroid_test_data")

    # We'll need an integrated reflection table and an experiment list.
    reflections_file = data_dir.join("integrated.pickle").strpath
    experiments_file = data_dir.join("experiments.json").strpath

    for arguments in (
        None,
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


def test_command_line(dials_data, tmpdir):
    """Test that the program runs correctly."""
    data_dir = dials_data("centroid_test_data")

    # We'll need an integrated reflection table and an experiment list.
    reflections_file = data_dir.join("integrated.pickle").strpath
    experiments_file = data_dir.join("experiments.json").strpath

    with tmpdir.as_cwd():
        run([experiments_file, reflections_file])

    output = tmpdir.join("corrected.refl")

    assert output.check(file=True)

    logfile = tmpdir.join("dials.anvil_correction.log").read()

    assert "Correcting integrated reflection intensities" in logfile
    assert "Writing the reflection table" in logfile
