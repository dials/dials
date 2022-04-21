"""
Tests for the dials.algorithms.scaling.plots module
"""

from __future__ import annotations

import numpy as np

from dials.algorithms.scaling.model.model import (
    ArrayScalingModel,
    PhysicalScalingModel,
    plot_scaling_models,
)
from dials.algorithms.scaling.plots import (
    normal_probability_plot,
    plot_array_absorption_plot,
    plot_array_decay_plot,
    plot_array_modulation_plot,
    plot_outliers,
)


def test_plot_array_absorption_plot():
    """Test plot generation for array absorption_correction"""
    array_dict = {
        "__id__": "array",
        "is_scaled": True,
        "absorption": {
            "n_parameters": 45,
            "parameters": list(range(1, 46)),
            "est_standard_devs": [],
        },
        "configuration_parameters": {
            "corrections": ["absorption"],
            "time_rot_interval": 10.0,
            "n_x_param": 3,
            "n_y_param": 3,
            "n_time_param": 5,
            "xmin": 0,
            "x_bin_width": 2,
            "ymin": 1,
            "y_bin_width": 2,
            "valid_osc_range": [0, 100],
        },
    }

    array_model = ArrayScalingModel.from_dict(array_dict)
    plot = plot_array_absorption_plot(array_model)
    assert plot["array_absorption_plot"]["data"][0]["x"]
    assert plot["array_absorption_plot"]["data"][0]["y"]
    assert plot["array_absorption_plot"]["data"][0]["z"]


def test_plot_array_decay_plot():
    """Test plot generation for array decay correction"""
    array_dict = {
        "__id__": "array",
        "is_scaled": True,
        "decay": {
            "n_parameters": 20,
            "parameters": list(range(1, 21)),
            "est_standard_devs": [],
        },
        "configuration_parameters": {
            "corrections": ["decay"],
            "time_rot_interval": 10.0,
            "n_res_param": 5,
            "res_bin_width": 0.1,
            "n_time_param": 4,
            "resmin": 0.05,
            "valid_osc_range": [0, 100],
        },
    }

    array_model = ArrayScalingModel.from_dict(array_dict)
    plot = plot_array_decay_plot(array_model)
    assert plot["array_decay_plot"]["data"][0]["x"]
    assert plot["array_decay_plot"]["data"][0]["y"]


def test_plot_array_modulation_plot():
    """Test plot generation for array modulation correction"""
    array_dict = {
        "__id__": "array",
        "is_scaled": True,
        "modulation": {
            "n_parameters": 25,
            "parameters": list(range(1, 26)),
            "est_standard_devs": [],
        },
        "configuration_parameters": {
            "corrections": ["modulation", "decay"],
            "n_x_mod_param": 5,
            "n_y_mod_param": 5,
            "xmin": 0,
            "x_det_bin_width": 2,
            "ymin": 1,
            "y_det_bin_width": 2,
            "time_rot_interval": 10.0,
            "n_res_param": 5,
            "res_bin_width": 0.1,
            "n_time_param": 4,
            "resmin": 0.05,
            "valid_osc_range": [0, 100],
        },
        "decay": {
            "n_parameters": 20,
            "parameters": list(range(1, 21)),
            "est_standard_devs": [],
        },
    }

    array_model = ArrayScalingModel.from_dict(array_dict)
    plot = plot_array_modulation_plot(array_model)
    assert plot["array_modulation_plot"]["data"][0]["x"]
    assert plot["array_modulation_plot"]["data"][0]["y"]


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
    d = plot_scaling_models(PhysicalScalingModel.from_dict(physical_dict))
    assert "smooth_scale_model" in d
    assert "absorption_surface" in d
    assert "absorption_parameters" in d
    assert d["smooth_scale_model"]["data"][0] != []
    assert d["absorption_parameters"]["data"][0] != []


def test_normal_probability_plot():
    data = {"delta_hl": np.arange(20, dtype=float)}
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
