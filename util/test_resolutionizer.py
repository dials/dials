from __future__ import absolute_import, division, print_function

import json
import pytest

from cctbx import uctbx
import iotbx.merging_statistics
from scitbx.array_family import flex
from scitbx.math import curve_fitting

from dials.command_line import resolutionizer as cmdline
from dials.util import resolution_analysis


def test_polynomial_fit():
    x = flex.double(range(-50, 50))
    p = (2, 3, 5)
    yo = flex.double(x.size())
    for i in range(len(p)):
        yo += p[i] * flex.pow(x, i)
    yf = resolution_analysis.polynomial_fit(x, yo, degree=2)
    assert yo == pytest.approx(yf)


def test_log_fit():
    x = flex.double(range(0, 100)) * 0.01
    p = (1, 2)
    yo = flex.double(x.size())
    for i in range(len(p)):
        yo += flex.exp(p[i] * flex.pow(x, i))
    yf = resolution_analysis.log_fit(x, yo, degree=2)
    assert yo == pytest.approx(yf, abs=1e-2)


def test_log_inv_fit():
    x = flex.double(range(0, 100)) * 0.01
    p = (1, 2)
    yo = flex.double(x.size())
    for i in range(len(p)):
        yo += 1 / flex.exp(p[i] * flex.pow(x, i))
    yf = resolution_analysis.log_inv_fit(x, yo, degree=2)
    assert yo == pytest.approx(yf, abs=1e-2)


def test_tanh_fit():
    x = flex.double(range(0, 100)) * 0.01
    f = curve_fitting.tanh(0.5, 1.5)
    yo = f(x)
    yf = resolution_analysis.tanh_fit(x, yo)
    assert yo == pytest.approx(yf, abs=1e-5)


@pytest.fixture
def merging_stats(dials_data):
    mtz = dials_data("x4wide_processed").join("AUTOMATIC_DEFAULT_scaled_unmerged.mtz")
    i_obs, _ = resolution_analysis.miller_array_from_mtz(mtz.strpath)
    return iotbx.merging_statistics.dataset_statistics(
        i_obs=i_obs,
        n_bins=20,
        binning_method="counting_sorted",
        use_internal_variance=False,
        eliminate_sys_absent=False,
        assert_is_not_unique_set_under_symmetry=False,
    )


def test_resolution_fit(merging_stats):
    d_star_sq = flex.double(uctbx.d_as_d_star_sq(b.d_min) for b in merging_stats.bins)
    y_obs = flex.double(b.r_merge for b in merging_stats.bins)
    result = resolution_analysis.resolution_fit(
        d_star_sq, y_obs, resolution_analysis.log_inv_fit, 0.6
    )
    assert result.d_min == pytest.approx(1.278, abs=1e-3)
    assert flex.max(flex.abs(result.y_obs - result.y_fit)) < 0.05


def test_resolution_cc_half(merging_stats):
    result = resolution_analysis.resolution_cc_half(merging_stats, limit=0.82)
    assert result.d_min == pytest.approx(1.242, abs=1e-3)
    result = resolution_analysis.resolution_cc_half(
        merging_stats,
        limit=0.82,
        cc_half_method="sigma_tau",
        model=resolution_analysis.polynomial_fit,
    )
    assert result.d_min == pytest.approx(1.23, abs=1e-3)
    assert flex.max(flex.abs(result.y_obs - result.y_fit)) < 0.04


def test_resolution_fit_from_merging_stats(merging_stats):
    result = resolution_analysis.resolution_fit_from_merging_stats(
        merging_stats, "i_over_sigma_mean", resolution_analysis.log_fit, limit=1.5
    )
    assert result.d_min == pytest.approx(1.295, abs=1e-3)
    assert flex.max(flex.abs(result.y_obs - result.y_fit)) < 1


def test_plot_result(merging_stats):
    result = resolution_analysis.resolution_cc_half(merging_stats, limit=0.82)
    d = resolution_analysis.plot_result("cc_half", result)
    assert "data" in d
    assert "layout" in d

    result = resolution_analysis.resolution_fit_from_merging_stats(
        merging_stats,
        "unmerged_i_over_sigma_mean",
        resolution_analysis.log_fit,
        limit=0.82,
    )
    d = resolution_analysis.plot_result("isigma", result)
    assert "data" in d
    assert "layout" in d


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
        "Resolution completeness:  1.20",
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


def test_resolutionizer_multi_sequence_with_batch_range(
    dials_data, run_in_tmpdir, capsys
):
    location = dials_data("l_cysteine_4_sweeps_scaled")
    refls = location.join("scaled_20_25.refl")
    expts = location.join("scaled_20_25.expt")

    cmdline.run(["batch_range=1900,3600", refls.strpath, expts.strpath],)
    captured = capsys.readouterr()

    expected_output = (
        "Resolution cc_half:       0.61",
        "Resolution I/sig:         0.59",
        "Resolution Mn(I/sig):     0.59",
    )
    for line in expected_output:
        assert line in captured.out
    assert run_in_tmpdir.join("dials.resolutionizer.html").check(file=1)
