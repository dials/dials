"""
Tests for the dials.algorithms.scaling.plots module
"""
from __future__ import absolute_import, division, print_function

from dials.algorithms.scaling.plots import (
    plot_scaling_models,
    plot_outliers,
    normal_probability_plot,
)


def test_plot_scaling_models():

    physical_dict = {
        "__id__": "physical",
        "is_scaled": True,
        "scale": {
            "n_parameters": 2,
            "parameters": [0.5, 1.0],
            "est_standard_devs": [0.05, 0.1],
        },
        "configuration_parameters": {
            "corrections": ["scale", "decay", "absorption"],
            "s_norm_fac": 0.1,
            "d_norm_fac": 0.1,
            "scale_rot_interval": 10.0,
            "decay_rot_interval": 10.0,
            "decay_restaint": 1e-1,
            "valid_osc_range": [0.0, 2.0],
        },
        "decay": {
            "n_parameters": 2,
            "parameters": [0.5, 1.0],
            "est_standard_devs": [0.05, 0.1],
        },
        "absorption": {
            "n_parameters": 4,
            "parameters": [0.1, -0.1, 0.05, -0.05],
            "est_standard_devs": [0.005, 0.005, 0.005, 0.005],
        },
    }
    d = plot_scaling_models(physical_dict)
    assert "smooth_scale_model" in d
    assert "absorption_surface" in d
    assert "absorption_parameters" in d
    assert d["smooth_scale_model"]["data"][0] != []
    assert d["absorption_parameters"]["data"][0] != []


def test_normal_probability_plot():
    data = {"delta_hl": list(range(20))}
    d = normal_probability_plot(data)
    assert "normal_distribution_plot" in d


def test_plot_outliers():
    """Test outlier plot, for standard and null data."""
    data = {"x": [1.0, 2.0], "y": [0.0, 1.0], "z": [1.0, 1.0], "image_size": [100, 200]}
    d = plot_outliers(data)
    assert "outliers_vs_z" in d
    assert "outlier_xy_positions" in d

    data = {"x": [], "y": [], "z": [], "image_size": [100, 200]}
    d = plot_outliers(data)
    assert "outliers_vs_z" in d
    assert "outlier_xy_positions" in d
