from __future__ import absolute_import, division, print_function

import os
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
    import libtbx.load_env

    paths = [
        os.path.join(libtbx.env.find_in_repositories("data-files/x4wide"), p)
        for p in input_files
    ]
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
Resolution I/sig:        1.53
Resolution Mn(I/sig):    1.51
Resolution Mn(I)/Mn(sig):    1.50"""
    for line in expected_output.splitlines():
        assert line in result.stdout

    expected_png = ["cc_half.png", "isigma.png", "misigma.png"]
    for png in expected_png:
        assert tmpdir.join(png).check()
