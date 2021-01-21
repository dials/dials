from unittest import mock

import pytest

from dxtbx.serialize import load

from dials.algorithms.scaling.observers import (
    ScalingHTMLContextManager,
    make_error_model_plots,
    make_filtering_plots,
    make_merging_stats_plots,
    make_outlier_plots,
    make_scaling_model_plots,
    print_scaling_model_error_summary,
    print_scaling_summary,
)
from dials.algorithms.scaling.scaling_library import scaled_data_as_miller_array
from dials.array_family import flex


@pytest.fixture
def test_data(dials_data):
    location = dials_data("l_cysteine_4_sweeps_scaled")
    refl = location.join("scaled_20_25.refl").strpath
    expt = location.join("scaled_20_25.expt").strpath

    refls = flex.reflection_table.from_file(refl).split_by_experiment_id()
    expts = load.experiment_list(expt, check_format=False)
    params = mock.Mock()
    return refls, expts, params


@pytest.fixture
def test_script(test_data):
    script = mock.Mock()
    refls, expts, params = test_data
    script.reflections = refls
    script.experiments = expts
    script.params = params
    script.scaled_miller_array = scaled_data_as_miller_array(refls, expts)
    script.merging_statistics_result = None
    script.anom_merging_statistics_result = None
    script.filtering_results = None
    return script


# test things which require refls, expts, params
def test_make_scaling_model_plots(test_data):
    _, expts, __ = test_data
    graphs = make_scaling_model_plots(expts)
    assert graphs
    assert len(graphs["scaling_model"]) == 6  # 2 models, 3 components each


def test_print_scaling_model_error_summary(test_data):
    _, expts, __ = test_data
    msg = print_scaling_model_error_summary(expts)
    assert (
        msg
        == """Warning: Over half (83.67%) of model parameters have signficant
uncertainty (sigma/abs(parameter) > 0.5), which could indicate a
poorly-determined scaling problem or overparameterisation.
"""
    )


def test_make_outlier_plots(test_data):
    refls, expts, _ = test_data
    graphs = make_outlier_plots(refls, expts)
    assert graphs
    assert len(graphs["outlier_plots"]) == 4


def test_make_error_model_plots(test_data):
    _, expts, params = test_data
    graphs = make_error_model_plots(params, expts)
    assert graphs
    assert "error_model_plots" in graphs


# test things which require a script object
def test_make_filtering_plots(test_script):
    plots = make_filtering_plots(test_script)
    assert plots


def test_make_merging_stats_plots(test_script):
    plots = make_merging_stats_plots(test_script)
    assert plots
    assert plots["image_range_tables"]


def test_print_scaling_summary(test_script):
    print_scaling_summary(test_script)


def test_ScalingHTMLContextManager(test_script, tmpdir):
    script = test_script
    script.params.output.html = tmpdir.join("test.html").strpath
    script.params.output.json = tmpdir.join("test.json").strpath
    with ScalingHTMLContextManager(script):
        pass

    assert tmpdir.join("test.html").check()
    assert tmpdir.join("test.json").check()
