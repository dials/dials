"""
Tests for the dials.report.plots module.
"""
import random
import pytest
from cctbx import miller, crystal
from scitbx.array_family import flex
from iotbx.merging_statistics import dataset_statistics
from dials.util.batch_handling import batch_manager
from dials.report.plots import (
    statistics_tables,
    cc_one_half_plot,
    i_over_sig_i_plot,
    i_over_sig_i_vs_batch_plot,
    scale_rmerge_vs_batch_plot,
)

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

def test_statistics_tables(iobs):
    """Test generation of statistics tables"""
    result = dataset_statistics(iobs, assert_is_not_unique_set_under_symmetry=False)
    tables = statistics_tables(result)
    assert len(tables) == 2  # overall and per resolution

def test_resolution_dependent_plots(iobs):
    """Test cc half plot, for centric and acentric data"""
    n_bins = 2
    result = dataset_statistics(
        iobs,
        assert_is_not_unique_set_under_symmetry=False,
        n_bins=n_bins
    )
    d = cc_one_half_plot(result)
    assert len(d['cc_one_half']['data']) == 4
    assert all([len(x['x']) == n_bins for x in d['cc_one_half']['data']])
    # check for correct labels
    assert list(d['cc_one_half']['layout']['xaxis']['ticktext']) == \
        ['1.26', '1.19', '1.13', '1.08', '1.04']

    assert list(d['cc_one_half']['layout']['xaxis']['tickvals']) == pytest.approx(
        [0.6319, 0.7055, 0.7792, 0.8528, 0.9264], 1e-4)

    # now try sigma tau and centric options
    d = cc_one_half_plot(result, method='sigma_tau', is_centric=True)
    assert len(d['cc_one_half']['data']) == 4
    assert all([len(x['x']) == n_bins for x in d['cc_one_half']['data'][:2]])
    assert d['cc_one_half']['data'][2] == {} # no anomalous plots
    assert d['cc_one_half']['data'][3] == {} # no anomalous plots
    # check for correct labels
    assert list(d['cc_one_half']['layout']['xaxis']['ticktext']) == \
        ['1.26', '1.19', '1.13', '1.08', '1.04']

    assert list(d['cc_one_half']['layout']['xaxis']['tickvals']) == pytest.approx(
        [0.6319, 0.7055, 0.7792, 0.8528, 0.9264], 1e-4)

    d = i_over_sig_i_plot(result)
    assert len(d['i_over_sig_i']['data'][0]['y']) == n_bins
    assert list(d['i_over_sig_i']['layout']['xaxis']['ticktext']) == \
        ['1.26', '1.19', '1.13', '1.08', '1.04']
    assert list(d['i_over_sig_i']['layout']['xaxis']['tickvals']) == pytest.approx(
        [0.6319, 0.7055, 0.7792, 0.8528, 0.9264], 1e-4)

@pytest.fixture
def batch_manager_fix():
    """Make a batch manager fixture"""

    batch_params = [{"id" : 0, "range" : [0, 10]}, {"id" : 1, "range" : [100, 110]}]
    batches = flex.int(range(0, 10) + range(100, 110))
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
