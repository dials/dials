from __future__ import annotations

import numpy as np


def flex_double_as_string(flex_array, n_digits=None):
    if n_digits is not None:
        flex_array = flex_array.round(n_digits)
    return list(flex_array.as_string())


def heatmap_unit_cell_scatter_plots(uc_params, nbins=100, mask_zeros=True):
    """
    Generate heatmap-style unit cell scatter plots.

    If mask_zeros is True, the output will have NaN values where there is no
    data, which will be absent (white) in the resulting plotly plots.
    """

    a, b, c = uc_params[0:3]

    def uc_param_hist2d(p1, p2):
        H, xedges, yedges = np.histogram2d(p1, p2, bins=nbins)
        H = np.rot90(H)
        H = np.flipud(H)
        if mask_zeros:
            Hmasked = np.ma.masked_where(H == 0, H)
            return xedges, yedges, Hmasked
        return xedges, yedges, H

    x1, y1, z1 = uc_param_hist2d(a, b)
    x2, y2, z2 = uc_param_hist2d(b, c)
    x3, y3, z3 = uc_param_hist2d(c, a)
    maxz = max(np.max(z1), np.max(z2), np.max(z3))
    return {
        "data": [
            {
                "x": x1.tolist(),
                "y": y1.tolist(),
                "z": z1.tolist(),
                "type": "heatmap",
                "mode": "markers",
                "name": "a vs. b",
                "colorscale": "Viridis",
                "xaxis": "x",
                "yaxis": "y",
                "zaxis": "z",
                "zmin": 1,
                "zmax": maxz,
            },
            {
                "x": x2.tolist(),
                "y": y2.tolist(),
                "z": z2.tolist(),
                "type": "heatmap",
                "mode": "markers",
                "name": "b vs. c",
                "colorscale": "Viridis",
                "xaxis": "x2",
                "yaxis": "y2",
                "zaxis": "z",
                "zmin": 1,
                "zmax": maxz,
            },
            {
                "x": x3.tolist(),
                "y": y3.tolist(),
                "z": z3.tolist(),
                "type": "heatmap",
                "mode": "markers",
                "name": "c vs. a",
                "colorbar": {
                    "title": "Frequency",
                    "titleside": "right",
                },
                "colorscale": "Viridis",
                "xaxis": "x3",
                "yaxis": "y3",
                "zmin": 1,
                "zmax": maxz,
            },
        ],
        "layout": {
            "title": "Distribution of unit cell parameters",
            "showlegend": False,
            "grid": {"rows": 1, "columns": 3, "subplots": [["xy", "x2y2", "x3y3"]]},
            "xaxis": {"title": "a (Å)"},
            "yaxis": {"title": "b (Å)"},
            "xaxis2": {"title": "b (Å)"},
            "yaxis2": {"title": "c (Å)"},
            "xaxis3": {"title": "c (Å)"},
            "yaxis3": {"title": "a (Å)"},
        },
        "help": """\
The distribution of the unit cell parameters: a vs b, b vs c and c vs a respectively.
""",
    }


def plot_uc_histograms(uc_params, scatter_style="points"):
    a, b, c, al, be, ga = (flex_double_as_string(p, n_digits=4) for p in uc_params)
    if scatter_style == "points":
        uc_scatter_plots = {
            "data": [
                {
                    "x": a,
                    "y": b,
                    "type": "scatter",
                    "mode": "markers",
                    "name": "a vs. b",
                    "xaxis": "x",
                    "yaxis": "y",
                },
                {
                    "x": b,
                    "y": c,
                    "type": "scatter",
                    "mode": "markers",
                    "name": "b vs. c",
                    "xaxis": "x2",
                    "yaxis": "y2",
                },
                {
                    "x": c,
                    "y": a,
                    "type": "scatter",
                    "mode": "markers",
                    "name": "c vs. a",
                    "xaxis": "x3",
                    "yaxis": "y3",
                },
            ],
            "layout": {
                "grid": {"rows": 1, "columns": 3, "pattern": "independent"},
                "title": "Distribution of unit cell parameters",
                "showlegend": False,
                "xaxis": {"title": "a (Å)"},
                "yaxis": {"title": "b (Å)"},
                "xaxis2": {"title": "b (Å)"},
                "yaxis2": {"title": "c (Å)"},
                "xaxis3": {"title": "c (Å)"},
                "yaxis3": {"title": "a (Å)"},
            },
            "help": """\
    The distribution of the unit cell parameters: a vs. b, b vs. c and c vs.a respectively.
    """,
        }
    elif scatter_style == "heatmap":
        uc_scatter_plots = heatmap_unit_cell_scatter_plots(uc_params)
    else:
        raise ValueError("Bad input for scatter_type parameter")

    plots = {
        "uc_scatter": uc_scatter_plots,
        "uc_hist": {
            "data": [
                {
                    "x": a,
                    "type": "histogram",
                    "connectgaps": False,
                    "name": "uc_hist_a",
                    "nbins": "auto",
                    "xaxis": "x",
                    "yaxis": "y",
                },
                {
                    "x": b,
                    "type": "histogram",
                    "connectgaps": False,
                    "name": "uc_hist_b",
                    "nbins": "auto",
                    "xaxis": "x2",
                    "yaxis": "y",
                },
                {
                    "x": c,
                    "type": "histogram",
                    "connectgaps": False,
                    "name": "uc_hist_c",
                    "nbins": "auto",
                    "xaxis": "x3",
                    "yaxis": "y",
                },
            ],
            "layout": {
                "grid": {"rows": 1, "columns": 3, "subplots": [["xy", "x2y", "x3y"]]},
                "title": "Histogram of unit cell parameters",
                "showlegend": False,
                "xaxis": {"title": "a (Å)"},
                "yaxis": {"title": "Frequency"},
                "xaxis2": {"title": "b (Å)"},
                "xaxis3": {"title": "c (Å)"},
            },
            "help": """\
    Histograms of unit cell parameters, a, b and c.
    """,
        },
    }
    uc_angle_hist = {
        "data": [],
        "help": "Histograms of unit cell angles.",
        "layout": {
            "title": "Histogram of unit cell angles",
            "showlegend": False,
            "yaxis": {"title": "Frequency"},
        },
    }
    subplots = []

    n = 0
    # only plot if not symmetry constrained, or only one crystal
    if len(set(al)) > 1 or len(al) == 1:
        n += 1
        uc_angle_hist["data"].append(
            {
                "x": al,
                "type": "histogram",
                "connectgaps": False,
                "name": "uc_hist_alpha",
                "nbins": "auto",
                "xaxis": "x" if n == 1 else f"x{n}",
                "yaxis": "y",
            }
        )
        subplots.append("xy" if n == 1 else f"x{n}y")
        uc_angle_hist["layout"]["xaxis" if n == 1 else f"xaxis{n}"] = {
            "title": "alpha (°)"
        }
    if len(set(be)) > 1 or len(be) == 1:
        n += 1
        uc_angle_hist["data"].append(
            {
                "x": be,
                "type": "histogram",
                "connectgaps": False,
                "name": "uc_hist_beta",
                "nbins": "auto",
                "xaxis": "x" if n == 1 else f"x{n}",
                "yaxis": "y",
            }
        )
        subplots.append("xy" if n == 1 else f"x{n}y")
        uc_angle_hist["layout"]["xaxis" if n == 1 else f"xaxis{n}"] = {
            "title": "beta (°)"
        }
    if len(set(ga)) > 1 or len(ga) == 1:
        n += 1
        uc_angle_hist["data"].append(
            {
                "x": ga,
                "type": "histogram",
                "connectgaps": False,
                "name": "uc_hist_gamma",
                "nbins": "auto",
                "xaxis": "x" if n == 1 else f"x{n}",
                "yaxis": "y",
            }
        )
        subplots.append("xy" if n == 1 else f"x{n}y")
        uc_angle_hist["layout"]["xaxis" if n == 1 else f"xaxis{n}"] = {
            "title": "gamma (°)"
        }
    if n:
        grid = {"rows": 1, "columns": 3, "subplots": [subplots]}
        uc_angle_hist["layout"]["grid"] = grid
        plots["uc_angle_hist"] = uc_angle_hist
    return plots


def scipy_dendrogram_to_plotly_json(ddict, title, xtitle=None, ytitle=None, help=None):
    colors = {
        "b": "rgb(31, 119, 180)",
        "g": "rgb(44, 160, 44)",
        "o": "rgb(255, 127, 14)",
        "r": "rgb(214, 39, 40)",
    }

    dcoord = ddict["dcoord"]
    icoord = ddict["icoord"]
    color_list = ddict["color_list"]
    ivl = ddict["ivl"]

    data = []
    xticktext = []
    xtickvals = []

    for k in range(len(dcoord)):
        x = icoord[k]
        y = dcoord[k]

        if y[0] == 0:
            xtickvals.append(x[0])
        if y[3] == 0:
            xtickvals.append(x[3])

        data.append(
            {
                "x": x,
                "y": y,
                "marker": {"color": colors.get(color_list[k])},
                "mode": "lines",
            }
        )

    xtickvals = sorted(xtickvals)
    xticktext = ivl
    d = {
        "data": data,
        "layout": {
            "barmode": "group",
            "legend": {"x": 100, "y": 0.5, "bordercolor": "transparent"},
            "margin": {"r": 10},
            "showlegend": False,
            "title": title,
            "xaxis": {
                "showline": False,
                "showgrid": False,
                "showticklabels": True,
                "tickangle": 300,
                "title": xtitle,
                "titlefont": {"color": "none"},
                "type": "linear",
                "ticktext": xticktext,
                "tickvals": xtickvals,
                "tickorientation": "vertical",
            },
            "yaxis": {
                "showline": False,
                "showgrid": False,
                "showticklabels": True,
                "tickangle": 0,
                "title": ytitle,
                "type": "linear",
            },
            "hovermode": "closest",
        },
    }
    if help:
        d["help"] = help
    return d
