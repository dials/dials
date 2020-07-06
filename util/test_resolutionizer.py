from __future__ import absolute_import, division, print_function

import pytest

from dials.command_line import resolutionizer


@pytest.mark.parametrize(
    "input_files",
    [
        ("AUTOMATIC_DEFAULT_scaled_unmerged.mtz",),
        ("AUTOMATIC_DEFAULT_scaled.refl", "AUTOMATIC_DEFAULT_scaled.expt"),
    ],
)
def test_resolutionizer(input_files, dials_data, run_in_tmpdir, capsys):
    paths = [dials_data("x4wide_processed").join(p).strpath for p in input_files]
    reference_mtz = dials_data("x4wide_processed").join("AUTOMATIC_DEFAULT_scaled.mtz")
    resolutionizer.run(
        [
            "plot=True",
            "cc_half=0.9",
            "isigma=2",
            "misigma=3",
            "rmerge=0.5",
            "completeness=1.0",
            "i_mean_over_sigma_mean=3",
            "batch_range=1,20",
            "batch_range=70,90",
            "space_group=P43212",
            "reference=%s" % reference_mtz,
            "cc_ref=0.9",
            "labels=IMEAN,SIGIMEAN",
        ]
        + paths,
    )
    captured = capsys.readouterr()
    expected_output = (
        "Resolution rmerge:        1.34",
        "Resolution completeness:  1.20",
        "Resolution cc_half:       1.56",
        "Resolution cc_ref:        1.3",
        "Resolution I/sig:         1.53",
        "Resolution Mn(I/sig):     1.51",
        "Resolution Mn(I)/Mn(sig): 1.50",
    )
    for line in expected_output:
        assert line in captured.out

    expected_png = (
        "cc_half.png",
        "isigma.png",
        "misigma.png",
        "completeness.png",
        "rmerge.png",
        "cc_ref.png",
        "i_mean_over_sigma_mean.png",
    )
    for png in expected_png:
        assert run_in_tmpdir.join(png).check()


def test_resolutionizer_multi_sequence_with_batch_range(
    dials_data, run_in_tmpdir, capsys
):
    location = dials_data("l_cysteine_4_sweeps_scaled")
    refls = location.join("scaled_20_25.refl")
    expts = location.join("scaled_20_25.expt")

    result = resolutionizer.run(
        ["batch_range=1900,3600", refls.strpath, expts.strpath],
    )
    captured = capsys.readouterr()
    assert not result.returncode and not result.stderr

    expected_output = (
        "Resolution cc_half:       0.61",
        "Resolution I/sig:         0.59",
        "Resolution Mn(I/sig):     0.59",
    )
    for line in expected_output:
        assert line in captured.out
