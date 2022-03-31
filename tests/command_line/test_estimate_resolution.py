from __future__ import annotations

import json

import procrunner
import pytest

from dials.command_line import estimate_resolution as cmdline


@pytest.mark.parametrize(
    "input_files",
    [
        ("AUTOMATIC_DEFAULT_scaled_unmerged.mtz",),
        ("AUTOMATIC_DEFAULT_scaled.refl", "AUTOMATIC_DEFAULT_scaled.expt"),
    ],
)
def test_x4wide(input_files, dials_data, run_in_tmp_path, capsys):
    x4wide = dials_data("x4wide_processed", pathlib=True)
    paths = [str(x4wide / p) for p in input_files]
    reference_mtz = x4wide / "AUTOMATIC_DEFAULT_scaled.mtz"
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
            f"reference={reference_mtz}",
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
    assert run_in_tmp_path.joinpath("resolutionizer.html").is_file()
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
    resolutionizer = run_in_tmp_path / "resolutionizer.json"
    assert resolutionizer.is_file()
    with resolutionizer.open() as fh:
        d = json.load(fh)
    assert set(d.keys()) == expected_keys


def test_multi_sequence_with_batch_range(dials_data, run_in_tmp_path, capsys):
    location = dials_data("l_cysteine_4_sweeps_scaled", pathlib=True)
    refls = location / "scaled_20_25.refl"
    expts = location / "scaled_20_25.expt"

    cmdline.run(
        ["batch_range=1900,3600", str(refls), str(expts)],
    )
    captured = capsys.readouterr()

    expected_output = "Resolution cc_half:       0.61"
    for line in expected_output:
        assert line in captured.out
    assert run_in_tmp_path.joinpath("dials.estimate_resolution.html").is_file()


def test_dispatcher_name():
    result = procrunner.run(["dials.estimate_resolution"])
    assert not result.returncode
    assert not result.stderr


def test_handle_fit_failure(dials_data, run_in_tmp_path, capsys):
    location = dials_data("l_cysteine_dials_output", pathlib=True)
    cmdline.run(
        ["misigma=1"]
        + [
            str(location / f)
            for f in (
                "11_integrated.expt",
                "11_integrated.refl",
                "23_integrated.expt",
                "23_integrated.refl",
            )
        ]
    )
    captured = capsys.readouterr()

    expected_output = (
        "Resolution fit against cc_half failed: No reflections left for fitting",
        "Resolution Mn(I/sig):     0.62",
    )
    for line in expected_output:
        assert line in captured.out
    assert run_in_tmp_path.joinpath("dials.estimate_resolution.html").is_file()


def test_mismatched_experiments_reflections(dials_data, run_in_tmp_path):
    location = dials_data("l_cysteine_dials_output", pathlib=True)
    with pytest.raises(SystemExit):
        cmdline.run(
            [
                str(location / f)
                for f in (
                    "11_integrated.expt",
                    "11_integrated.refl",
                    "23_integrated.refl",
                )
            ]
        )
