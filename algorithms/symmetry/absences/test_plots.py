"""Test the plotting module from systematic_absences."""
from dials.algorithms.symmetry.absences.plots import plot_screw_axes, color_axis_data


def test_color_axis_data():
    """Test the helper function used to color the data."""
    colors = color_axis_data("41c", list(range(1, 9)))
    assert colors == [0, 0, 0, 1, 0, 0, 0, 1]
    colors = color_axis_data("42c", list(range(1, 9)))
    assert colors == [0, 1, 0, 1, 0, 1, 0, 1]
    colors = color_axis_data("31c", list(range(1, 7)))
    assert colors == [0, 0, 1, 0, 0, 1]
    colors = color_axis_data("21c", list(range(1, 7)))
    assert colors == [0, 1, 0, 1, 0, 1]
    colors = color_axis_data("61c", list(range(1, 7)))
    assert colors == [0, 0, 0, 0, 0, 1]
    colors = color_axis_data("62c", list(range(1, 7)))
    assert colors == [0, 0, 1, 0, 0, 1]
    colors = color_axis_data("63c", list(range(1, 7)))
    assert colors == [0, 1, 0, 1, 0, 1]


def test_plot_screw_axes():
    """Test that the plotting function returns plots with data."""
    data = {
        "miller_axis_vals": [1, 2, 3, 4, 5, 6],
        "i_over_sigma": [0.5, 20.0, 0.2, 10.0, 0.4, 15.0],
        "intensities": [0.5, 20.0, 0.2, 10.0, 0.4, 15.0],
        "sigmas": [1.0] * 6,
    }
    plots = plot_screw_axes({"21a": data})
    for pl in plots.values():
        assert pl["data"][0]["x"]
        assert pl["data"][0]["y"]
