"""Definitions of plots for systematic absences."""

from __future__ import annotations


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
    d = {}
    for name, data in screw_axes_data.items():
        d.update(
            {
                "plot_" + name: {
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
                        "title": f"I (merged) / σ (merged) along axis {name}",
                        "xaxis": {
                            "domain": [0, 1],
                            "anchor": "y",
                            "title": "index along axis",
                        },
                        "yaxis": {"domain": [0, 1], "anchor": "x", "title": "I/σ(I)"},
                    },
                }
            }
        )
        if data["fourier_space_data"]:
            # make frequency plots
            # need to add vertical lines to indicate axis repeats
            y = data["fourier_space_data"]["fourier_space"]
            y_min = min(y)
            y_max = max(y)
            n = len(y)
            xtickvals = [float(n // 2)]
            xticktext = ["1/2", "1/3", "2/3", "1/4", "3/4", "1/6", "5/6"]
            for i in [3, 4, 6]:
                xloc = float(n // i)
                xloc2 = n - xloc
                xtickvals.extend([xloc, xloc2])

            plot = {
                "frequencies_plot_" + name: {
                    "data": [
                        {
                            "x": list(range(len(y))),
                            "y": list(y),
                            "type": "scatter",
                            "name": "Fourier amplitudes",
                            "xaxis": "x",
                            "yaxis": "y",
                            "mode": "lines",
                        },
                    ],
                    "layout": {
                        "title": "Fourier amplitudes for axis " + name,
                        "xaxis": {
                            "title": "frequency",
                            "tickvals": xtickvals,
                            "ticktext": xticktext,
                        },
                        "yaxis": {"title": "Amplitude", "domain": [y_min, y_max]},
                    },
                    "help": """
This plot shows the amplitudes of the frequency spectrum of I/sigma, calculated via
a discrete fourier transform. The frequencies corresponding to systematic absences
are indicated by vertical lines.""",
                }
            }
            for i, dash in zip([2, 3, 4, 6], ["dash", "dot", "dashdot", "solid"]):
                xloc = float(n // i)
                plot["frequencies_plot_" + name]["data"].append(
                    {
                        "x": [xloc, xloc],
                        "y": [y_min, y_max],
                        "type": "scatter",
                        "xaxis": "x",
                        "yaxis": "y",
                        "mode": "lines",
                        "line": {
                            "color": "black" if i == data["axis_repeat"] else "grey",
                            "dash": dash,
                        },
                        "name": f"{i}-fold repeat",
                    }
                )
                if i == 2:
                    continue  # don't plot same thing twice
                xloc2 = float(n - xloc)
                plot["frequencies_plot_" + name]["data"].append(
                    {
                        "x": [xloc2, xloc2],
                        "y": [y_min, y_max],
                        "type": "scatter",
                        "xaxis": "x",
                        "yaxis": "y",
                        "mode": "lines",
                        "line": {
                            "color": "black" if i == data["axis_repeat"] else "grey",
                            "dash": dash,
                        },
                        "name": None,
                        "showlegend": False,
                    }
                )
            d.update(plot)

        d.update(
            {
                "intensities_plot_" + name: {
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
                        "title": f"I, σ (merged) along axis {name}",
                        "xaxis": {
                            "domain": [0, 1],
                            "anchor": "y",
                            "title": "index along axis",
                        },
                        "yaxis": {"domain": [0, 1], "anchor": "x", "title": "I, σ"},
                    },
                }
            }
        )
    return d
