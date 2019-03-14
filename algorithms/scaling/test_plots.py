from scitbx.array_family import flex
from iotbx.merging_statistics import dataset_statistics
from dials.algorithms.scaling.plots import (
    plot_scaling_models,
    plot_outliers,
    statistics_tables,
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


def test_normal_probability_plot():
    data = {"delta_hl": range(20)}
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


def test_statistics_tables():
    """Test generation of statistics tables"""
    import random
    from cctbx import miller
    from cctbx import crystal

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
    result = dataset_statistics(iobs, assert_is_not_unique_set_under_symmetry=False)
    tables = statistics_tables(result)
    assert len(tables) == 2  # overall and per resolution
