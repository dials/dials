"""Definitions of plots for systematic absences."""
from __future__ import absolute_import, division, print_function
from collections import OrderedDict


def color_axis_data(name, miller_axis_vals):
    """Generate a sequence for coloring the datapoints for a screw axis."""
    if name.startswith("41"):
        colors = [1 if m % 4 == 0 else 0 for m in miller_axis_vals]
    elif name.startswith("21") or name.startswith("42") or name.startswith("63"):
        colors = [1 if m % 2 == 0 else 0 for m in miller_axis_vals]
    elif name.startswith("31") or name.startswith("62"):
        colors = [1 if m % 3 == 0 else 0 for m in miller_axis_vals]
    elif name.startswith("61"):
        colors = [1 if m % 6 == 0 else 0 for m in miller_axis_vals]
    return colors


def plot_screw_axes(screw_axes_data):
    """Generate scatter plot data for screw axes."""
    d = OrderedDict()
    for name, data in screw_axes_data.items():
        d.update(
            {
                "plot_"
                + name: {
                    "data": [
                        {
                            "x": list(data["miller_axis_vals"]),
                            "y": list(data["i_over_sigma"]),
                            "type": "scatter",
                            "name": name,
                            "xaxis": "x",
                            "yaxis": "y",
                            "mode": "markers",
                            "marker": {
                                "color": color_axis_data(
                                    name, list(data["miller_axis_vals"])
                                ),
                                "colorscale": "Viridis",
                            },
                        }
                    ],
                    "layout": {
                        "title": "I (merged) / sigma (merged) along axis %s" % name,
                        "xaxis": {
                            "domain": [0, 1],
                            "anchor": "y",
                            "title": "index along axis",
                        },
                        "yaxis": {
                            "domain": [0, 1],
                            "anchor": "x",
                            "title": "I / sigma",
                        },
                    },
                }
            }
        )
        d.update(
            {
                "intensities_plot_"
                + name: {
                    "data": [
                        {
                            "x": list(data["miller_axis_vals"]),
                            "y": list(data["intensities"]),
                            "type": "scatter",
                            "name": "intensity",
                            "xaxis": "x",
                            "yaxis": "y",
                            "mode": "markers",
                            "marker": {
                                "color": color_axis_data(
                                    name, list(data["miller_axis_vals"])
                                ),
                                "colorscale": "Viridis",
                            },
                        },
                        {
                            "x": list(data["miller_axis_vals"]),
                            "y": list(data["sigmas"]),
                            "type": "scatter",
                            "name": "sigma",
                            "xaxis": "x",
                            "yaxis": "y",
                            "mode": "markers",
                        },
                    ],
                    "layout": {
                        "title": "I, sigma (merged) along axis %s" % name,
                        "xaxis": {
                            "domain": [0, 1],
                            "anchor": "y",
                            "title": "index along axis",
                        },
                        "yaxis": {"domain": [0, 1], "anchor": "x", "title": "I, sigma"},
                    },
                }
            }
        )
    return d
