import json

import pytest

from dials.command_line import estimate_resolution as cmdline


@pytest.mark.parametrize(
    "input_files",
    [
        ("AUTOMATIC_DEFAULT_scaled_unmerged.mtz",),
        ("AUTOMATIC_DEFAULT_scaled.refl", "AUTOMATIC_DEFAULT_scaled.expt"),
    ],
)
def test_x4wide(input_files, dials_data, run_in_tmpdir, capsys):
    paths = [dials_data("x4wide_processed").join(p).strpath for p in input_files]
    reference_mtz = dials_data("x4wide_processed").join("AUTOMATIC_DEFAULT_scaled.mtz")
    result = cmdline.run(
        [
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
            "html=resolutionizer.html",
            "json=resolutionizer.json",
        ]
        + paths,
    )
    captured = capsys.readouterr()
    expected_output = (
        "Resolution rmerge:        1.34",
        "Resolution cc_half:       1.56",
        "Resolution cc_ref:        1.3",
        "Resolution I/sig:         1.53",
        "Resolution Mn(I/sig):     1.51",
        "Resolution Mn(I)/Mn(sig): 1.50",
    )
    for line in expected_output:
        assert line in captured.out
    assert run_in_tmpdir.join("resolutionizer.html").check(file=1)
    expected_keys = {
        "cc_half",
        "cc_ref",
        "isigma",
        "misigma",
        "i_mean_over_sigma_mean",
        "rmerge",
        "completeness",
    }
    assert set(result.keys()) == expected_keys
    assert run_in_tmpdir.join("resolutionizer.json").check(file=1)
    with run_in_tmpdir.join("resolutionizer.json").open("r") as fh:
        d = json.load(fh)
    assert set(d.keys()) == expected_keys


def test_multi_sequence_with_batch_range(dials_data, run_in_tmpdir, capsys):
    location = dials_data("l_cysteine_4_sweeps_scaled")
    refls = location.join("scaled_20_25.refl")
    expts = location.join("scaled_20_25.expt")

    cmdline.run(
        ["batch_range=1900,3600", refls.strpath, expts.strpath],
    )
    captured = capsys.readouterr()

    expected_output = (
        "Resolution cc_half:       0.61",
        "Resolution I/sig:         0.59",
        "Resolution Mn(I/sig):     0.59",
    )
    for line in expected_output:
        assert line in captured.out
    assert run_in_tmpdir.join("dials.estimate_resolution.html").check(file=1)


@pytest.mark.xfail("os.name == 'nt'", reason="warnings do not go to stderr")
def test_dispatcher_name():
    import procrunner

    result = procrunner.run(["dials.resolutionizer"])
    assert not result.returncode
    assert (
        b"dials.resolutionizer is now deprecated, please use dials.estimate_resolution instead"
        in result.stderr
    )

    result = procrunner.run(["dials.estimate_resolution"])
    assert not result.returncode
    assert not result.stderr


def test_handle_fit_failure(dials_data, run_in_tmpdir, capsys):
    location = dials_data("l_cysteine_dials_output")
    filenames = [
        location.join("11_integrated.expt"),
        location.join("11_integrated.refl"),
        location.join("23_integrated.expt"),
        location.join("23_integrated.refl"),
    ]
    cmdline.run(["misigma=1"] + [f.strpath for f in filenames])
    captured = capsys.readouterr()

    expected_output = (
        "Resolution fit against cc_half failed: No reflections left for fitting",
        "Resolution Mn(I/sig):     0.62",
    )
    for line in expected_output:
        assert line in captured.out
    assert run_in_tmpdir.join("dials.estimate_resolution.html").check(file=1)
