# -*- coding: utf-8 -*-

from __future__ import absolute_import, division, print_function
from collections import OrderedDict


def flex_double_as_string(flex_array, n_digits=None):
    if n_digits is not None:
        flex_array = flex_array.round(n_digits)
    return list(flex_array.as_string())


def plot_uc_histograms(uc_params):

    a, b, c = (flex_double_as_string(p, n_digits=4) for p in uc_params[:3])
    d = OrderedDict()

    d["uc_scatter"] = {
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
    }

    d["uc_hist"] = {
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
    }

    return d


def scipy_dendrogram_to_plotly_json(ddict, title, xtitle=None, ytitle=None):
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
    return d
