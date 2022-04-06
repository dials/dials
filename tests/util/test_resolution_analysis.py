from __future__ import annotations

import pytest

import iotbx.merging_statistics
from cctbx import uctbx
from scitbx.array_family import flex
from scitbx.math import curve_fitting

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
    mtz = str(
        dials_data("x4wide_processed", pathlib=True)
        / "AUTOMATIC_DEFAULT_scaled_unmerged.mtz"
    )
    i_obs, _ = resolution_analysis.miller_array_from_mtz(mtz)
    return iotbx.merging_statistics.dataset_statistics(
        i_obs=i_obs,
        n_bins=20,
        binning_method="counting_sorted",
        use_internal_variance=False,
        eliminate_sys_absent=False,
        assert_is_not_unique_set_under_symmetry=False,
        cc_one_half_significance_level=0.1,
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
    assert result.d_min == pytest.approx(1.233, abs=1e-3)
    assert flex.max(flex.abs(result.y_obs - result.y_fit)) < 0.04
    assert result.critical_values is not None
    assert len(result.critical_values) == len(result.d_star_sq)


def test_resolution_fit_from_merging_stats(merging_stats):
    result = resolution_analysis.resolution_fit_from_merging_stats(
        merging_stats, "i_over_sigma_mean", resolution_analysis.log_fit, limit=1.5
    )
    assert result.d_min == pytest.approx(1.295, abs=1e-3)
    assert flex.max(flex.abs(result.y_obs - result.y_fit)) < 1


def test_resolution_fit_interpolation_error(merging_stats):
    result = resolution_analysis.resolution_fit_from_merging_stats(
        merging_stats, "i_over_sigma_mean", resolution_analysis.log_fit, limit=25
    )
    assert result.d_min is None


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
