from __future__ import absolute_import, division, print_function

import pytest

import procrunner


@pytest.mark.parametrize(
    "input_files",
    [
        ("AUTOMATIC_DEFAULT_scaled_unmerged.mtz",),
        ("AUTOMATIC_DEFAULT_scaled.refl", "AUTOMATIC_DEFAULT_scaled.expt"),
    ],
)
def test_resolutionizer(input_files, dials_data, tmpdir):
    paths = [dials_data("x4wide_processed").join(p) for p in input_files]
    reference_mtz = dials_data("x4wide_processed").join("AUTOMATIC_DEFAULT_scaled.mtz")
    result = procrunner.run(
        [
            "dials.resolutionizer",
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
        working_directory=tmpdir,
    )
    assert not result.returncode and not result.stderr
    print(result.stdout)

    expected_output = """\
Resolution rmerge:       1.34
Resolution completeness: 1.20
Resolution cc_half:      1.62
Resolution cc_ref:       1.31
Resolution I/sig:        1.53
Resolution Mn(I/sig):    1.51
Resolution Mn(I)/Mn(sig):    1.50"""
    for line in expected_output.encode("latin-1").splitlines():
        assert line in result.stdout

    expected_png = [
        "cc_half.png",
        "isigma.png",
        "misigma.png",
        "completeness.png",
        "rmerge.png",
        "cc_ref.png",
        "i_mean_over_sigma_mean.png",
    ]
    for png in expected_png:
        assert tmpdir.join(png).check()
