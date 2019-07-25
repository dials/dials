"""
Tests for the dials.report.plots module.
"""
from __future__ import absolute_import, division, print_function

import itertools
import random

import mock as mock
import pytest
from cctbx import miller, crystal
from cctbx.array_family import flex
from dials.report.plots import (
    ResolutionPlotsAndStats,
    i_over_sig_i_vs_batch_plot,
    scale_rmerge_vs_batch_plot,
    IntensityStatisticsPlots,
    i_over_sig_i_vs_i_plot,
    AnomalousPlotter,
)
from dials.util.batch_handling import batch_manager
from iotbx.merging_statistics import dataset_statistics


@pytest.fixture
def iobs():
    """Generate a miller array for testing plotters."""
    ms = miller.build_set(
        crystal_symmetry=crystal.symmetry(
            space_group_symbol="P1", unit_cell=(6, 6, 6, 90, 90, 90)
        ),
        anomalous_flag=True,
        d_min=1.0,
    )
    data = flex.double(float(random.randrange(0, 100)) for _ in range(ms.size()))
    iobs = miller.array(ms, data, sigmas=data)
    iobs.change_symmetry(space_group_symbol="P222", merge_non_unique=False)
    iobs.set_info(miller.array_info(source="DIALS", source_type="reflection_tables"))
    iobs.set_observation_type_xray_intensity()
    return iobs


def test_AnomalousPlotter():

    "Make a larger array to allow all plots to be made"
    cs = crystal.symmetry(space_group_symbol="P1", unit_cell=(6, 6, 6, 90, 90, 90))
    ms = miller.build_set(cs, anomalous_flag=True, d_min=1.0)
    indices = ms.indices()
    new_indices = flex.miller_index(list(indices) * 10)
    new_ms = miller.set(crystal_symmetry=cs, indices=new_indices, anomalous_flag=True)
    data = flex.double(float(random.randrange(1, 100)) for _ in range(new_ms.size()))
    iobs = miller.array(new_ms, data, sigmas=data)
    iobs.change_symmetry(space_group_symbol="P222", merge_non_unique=False)
    iobs.set_info(miller.array_info(source="DIALS", source_type="reflection_tables"))
    iobs.set_observation_type_xray_intensity()

    plotter = AnomalousPlotter(iobs)
    d = plotter.make_plots()
    expected = ["normal_distribution_plot_highres", "anom_correl_plot"]
    keys = d.keys()
    for k in expected:
        assert k in keys
        assert d[k]["data"][0]["x"]  # check some data there

    plotter = AnomalousPlotter(iobs, strong_cutoff=3.0)
    d = plotter.make_plots()
    expected = [
        "anom_correl_plot",
        "anom_scatter_plot_lowres",
        "normal_distribution_plot_lowres",
    ]
    keys = list(d.keys())
    for k in expected:
        assert k in keys
        assert d[k]["data"][0]["x"]  # check some data there


def test_IntensityStatisticsPlots(iobs):
    n_bins = 2
    # mock the wilson scaling
    wilson_scaling = mock.Mock()
    wilson_scaling.d_star_sq = [1.0, 2.0]
    wilson_scaling.mean_I_obs_data = [1.0, 2.0]
    wilson_scaling.mean_I_obs_theory = [1.0, 2.0]
    wilson_scaling.mean_I_normalisation = [1.0, 2.0]
    # mock the twin results
    twin_results = mock.Mock()
    twin_results.nz_test.z = [1.0, 2.0]
    twin_results.nz_test.ac_obs = [1.0, 2.0]
    twin_results.nz_test.c_obs = [1.0, 2.0]
    twin_results.nz_test.ac_untwinned = [1.0, 2.0]
    twin_results.nz_test.c_untwinned = [1.0, 2.0]
    twin_results.l_test.l_values = [1.0, 2.0]
    twin_results.l_test.l_cumul_untwinned = [1.0, 2.0]
    twin_results.l_test.l_cumul_perfect_twin = [1.0, 2.0]
    twin_results.l_test.l_cumul = [1.0, 2.0]
    # mock the xtraige analysis
    xtriage_analyses = mock.Mock()
    xtriage_analyses.wilson_scaling = wilson_scaling
    xtriage_analyses.twin_results = twin_results

    plotter = IntensityStatisticsPlots(
        iobs, n_resolution_bins=n_bins, xtriage_analyses=xtriage_analyses
    )

    d = plotter.generate_resolution_dependent_plots()
    assert "wilson_intensity_plot" in d
    assert "second_moments" in d
    d = plotter.generate_miscellanous_plots()
    assert "cumulative_intensity_distribution" in d
    assert "l_test" in d

    def mock_xtriage(*args, **kwargs):
        return xtriage_analyses

    with mock.patch("mmtbx.scaling.xtriage.xtriage_analyses", new=mock_xtriage):
        plotter = IntensityStatisticsPlots(iobs)
        d = plotter.generate_resolution_dependent_plots()
        assert "wilson_intensity_plot" in d
        assert "second_moments" in d
        d = plotter.generate_miscellanous_plots()
        assert "cumulative_intensity_distribution" in d
        assert "l_test" in d

    # try with anomalous
    plotter = IntensityStatisticsPlots(
        iobs,
        anomalous=True,
        n_resolution_bins=n_bins,
        xtriage_analyses=xtriage_analyses,
    )

    d = plotter.generate_resolution_dependent_plots()
    assert "wilson_intensity_plot" in d
    assert "second_moments" in d
    d = plotter.generate_miscellanous_plots()
    assert "cumulative_intensity_distribution" in d
    assert "l_test" in d
    assert "multiplicities" in d


def test_i_over_sig_i_vs_i_plot(iobs):
    """Test the generation of 2d heatmap plots for intensity dist."""
    d = i_over_sig_i_vs_i_plot(iobs.data(), iobs.sigmas())
    assert "i_over_sig_isq_vs_i" in d
    assert "i_over_sig_isq_vs_i" in d


def test_ResolutionPlotsAndStats(iobs):
    i_obs_anom = iobs.as_anomalous_array()
    iobs_anom = i_obs_anom.map_to_asu().customized_copy(info=iobs.info())
    n_bins = 2
    result = dataset_statistics(
        iobs, assert_is_not_unique_set_under_symmetry=False, n_bins=n_bins
    )
    anom_result = dataset_statistics(
        iobs_anom,
        assert_is_not_unique_set_under_symmetry=False,
        anomalous=True,
        n_bins=n_bins,
    )
    plotter = ResolutionPlotsAndStats(result, anom_result)

    assert plotter.d_star_sq_ticktext == ["1.26", "1.19", "1.13", "1.08", "1.04"]

    assert plotter.d_star_sq_tickvals == pytest.approx(
        [0.6319, 0.7055, 0.7792, 0.8528, 0.9264], 1e-4
    )

    tables = plotter.statistics_tables()
    assert len(tables) == 2  # overall and per resolution

    # test plots individually
    d = plotter.cc_one_half_plot()
    assert len(d["cc_one_half"]["data"]) == 4
    assert all([len(x["x"]) == n_bins for x in d["cc_one_half"]["data"]])

    d = plotter.i_over_sig_i_plot()
    assert len(d["i_over_sig_i"]["data"]) == 1
    assert len(d["i_over_sig_i"]["data"][0]["y"]) == n_bins

    d = plotter.completeness_plot()
    assert len(d["completeness"]["data"]) == 2
    assert len(d["completeness"]["data"][0]["y"]) == n_bins

    d = plotter.multiplicity_vs_resolution_plot()
    assert len(d["multiplicity_vs_resolution"]["data"]) == 2
    assert len(d["multiplicity_vs_resolution"]["data"][0]["y"]) == n_bins

    # now try centric options and sigma tau for cc_one_half
    plotter = ResolutionPlotsAndStats(result, anom_result, is_centric=True)
    d = plotter.cc_one_half_plot(method="sigma_tau")
    assert len(d["cc_one_half"]["data"]) == 4
    assert all([len(x["x"]) == n_bins for x in d["cc_one_half"]["data"][:2]])
    assert d["cc_one_half"]["data"][2] == {}  # no anomalous plots
    assert d["cc_one_half"]["data"][3] == {}  # no anomalous plots
    d = plotter.completeness_plot()
    assert len(d["completeness"]["data"]) == 2
    assert len(d["completeness"]["data"][0]["y"]) == n_bins
    assert d["completeness"]["data"][1] == {}
    d = plotter.multiplicity_vs_resolution_plot()
    assert len(d["multiplicity_vs_resolution"]["data"]) == 2
    assert len(d["multiplicity_vs_resolution"]["data"][0]["y"]) == n_bins
    assert d["multiplicity_vs_resolution"]["data"][1] == {}

    plots = plotter.make_all_plots()
    for plot in plots.values():
        assert plot["layout"]["xaxis"]["ticktext"] == plotter.d_star_sq_ticktext
        assert plot["layout"]["xaxis"]["tickvals"] == plotter.d_star_sq_tickvals


@pytest.fixture
def batch_manager_fix():
    """Make a batch manager fixture"""

    batch_params = [{"id": 0, "range": [0, 10]}, {"id": 1, "range": [100, 110]}]
    batches = flex.int(itertools.chain(range(0, 10), range(100, 110)))
    return batch_manager(batches, batch_params)


def test_i_over_sig_i_vs_batch_plot(batch_manager_fix):
    """Test the IsigI batch plot"""
    bm = batch_manager_fix
    isigi = flex.double(range(0, 20))
    d = i_over_sig_i_vs_batch_plot(bm, isigi)
    assert list(d["i_over_sig_i_vs_batch"]["data"][0]["x"]) == list(bm.reduced_batches)
    assert list(d["i_over_sig_i_vs_batch"]["data"][0]["y"]) == list(isigi)


def test_scale_rmerge_vs_batch_plot(batch_manager_fix):
    """Test the scale and rmerge batch plot. Should have the option
    to plot without scales (for xia2)."""
    bm = batch_manager_fix
    rmergevsb = flex.double(range(0, 20))
    scalesvsb = flex.double(range(1, 21))
    d = scale_rmerge_vs_batch_plot(bm, rmergevsb, scalesvsb)
    assert list(d["scale_rmerge_vs_batch"]["data"][0]["x"]) == list(bm.reduced_batches)
    assert list(d["scale_rmerge_vs_batch"]["data"][1]["x"]) == list(bm.reduced_batches)
    assert list(d["scale_rmerge_vs_batch"]["data"][0]["y"]) == list(scalesvsb)
    assert list(d["scale_rmerge_vs_batch"]["data"][1]["y"]) == list(rmergevsb)

    d = scale_rmerge_vs_batch_plot(bm, rmergevsb)
    assert d["scale_rmerge_vs_batch"]["data"][0] == {}
    assert list(d["scale_rmerge_vs_batch"]["data"][1]["x"]) == list(bm.reduced_batches)
    assert d["scale_rmerge_vs_batch"]["data"][0] == {}
    assert list(d["scale_rmerge_vs_batch"]["data"][1]["y"]) == list(rmergevsb)
