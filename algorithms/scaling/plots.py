# -*- coding: utf-8 -*-
"""
Make plotly plots for html output by dials.scale, dials.report or xia2.report.
"""
from __future__ import absolute_import, division, print_function
from collections import OrderedDict
import math as pymath
import numpy as np
from scitbx import math as scitbxmath
from scitbx.math import distributions
from dials.array_family import flex
from dials.algorithms.scaling.model.model import PhysicalScalingModel


def plot_scaling_models(scaling_model_dict):
    d = OrderedDict()
    if scaling_model_dict["__id__"] == "physical":
        model = PhysicalScalingModel.from_dict(scaling_model_dict)
        d.update(_plot_smooth_scales(model))
        if "absorption" in model.components:
            d.update(plot_absorption_parameters(model))
            d.update(plot_absorption_surface(model))
    return d


def _get_smooth_plotting_data_from_model(physical_model, component="scale"):
    """Return a tuple of phis, parameters, parameter esds,
    sample positions for plotting and sample scale values."""
    configdict = physical_model.configdict
    valid_osc = configdict["valid_osc_range"]
    sample_values = flex.double(
        np.linspace(
            valid_osc[0],
            valid_osc[1],
            int((valid_osc[1] - valid_osc[0]) / 0.1) + 1,
            endpoint=True,
        )
    )  # Make a grid of
    # points with 10 points per degree.
    if component == "scale":
        if "scale" in configdict["corrections"]:
            scale_SF = physical_model.components["scale"]
            norm_vals = sample_values / physical_model.configdict["scale_rot_interval"]
            scale_SF.data = {"x": norm_vals}
            scale_SF.update_reflection_data()
            s = scale_SF.calculate_scales()
            smoother_phis = [
                (i * configdict["scale_rot_interval"]) + valid_osc[0]
                for i in scale_SF.smoother.positions()
            ]
        return (
            smoother_phis,
            scale_SF.parameters,
            scale_SF.parameter_esds,
            sample_values,
            s,
        )
    if component == "decay":
        if "decay" in configdict["corrections"]:
            decay_SF = physical_model.components["decay"]
            norm_vals = sample_values / physical_model.configdict["decay_rot_interval"]
            decay_SF.data = {"x": norm_vals, "d": flex.double(norm_vals.size(), 1.0)}
            decay_SF.update_reflection_data()
            s = decay_SF.calculate_scales()
            smoother_phis = [
                (i * configdict["decay_rot_interval"]) + valid_osc[0]
                for i in decay_SF._smoother.positions()
            ]
        return (
            smoother_phis,
            decay_SF.parameters,
            decay_SF.parameter_esds,
            sample_values,
            s,
        )


smooth_help_msg = """
The inverse scale factor g, for a given reflection i, is defined as the product
of the individual model components.
For the physical scaling model, up to three terms are defined:
smoothly-varying scaling term Ci,
isotropic radiation damage term Ri = exp( Bi / (2 * di * di) ),
absorption surface correction Si,
The scaled intensity is therefore given by Ii / gi, where
gi = Ci * Ri * Si
"""


def _plot_smooth_scales(physical_model):
    """Plot smooth scale factors for the physical model."""

    d = {
        "smooth_scale_model": {
            "data": [],
            "layout": {
                "title": "Smoothly varying corrections",
                "xaxis": {
                    "domain": [0, 1],
                    "anchor": "y",
                    "title": "rotation angle (degrees)",
                },
                "yaxis": {
                    "domain": [0, 0.45],
                    "anchor": "x",
                    "title": "relative B-factor",
                },
                "yaxis2": {
                    "domain": [0.5, 0.95],
                    "anchor": "x",
                    "title": "inverse <br> scale factor",
                },
            },
            "help": smooth_help_msg,
        }
    }

    data = []

    if "scale" in physical_model.components:
        smoother_phis, parameters, parameter_esds, sample_values, sample_scales = _get_smooth_plotting_data_from_model(
            physical_model, component="scale"
        )

        data.append(
            {
                "x": list(sample_values),
                "y": list(sample_scales),
                "type": "line",
                "name": "smoothly-varying <br>scale correction",
                "xaxis": "x",
                "yaxis": "y2",
            }
        )
        data.append(
            {
                "x": list(smoother_phis),
                "y": list(parameters),
                "type": "scatter",
                "mode": "markers",
                "name": "smoothly-varying <br>scale parameters",
                "xaxis": "x",
                "yaxis": "y2",
            }
        )
        if parameter_esds:
            data[-1]["error_y"] = {"type": "data", "array": list(parameter_esds)}

    if "decay" in physical_model.components:
        smoother_phis, parameters, parameter_esds, sample_values, sample_scales = _get_smooth_plotting_data_from_model(
            physical_model, component="decay"
        )

        data.append(
            {
                "x": list(sample_values),
                "y": list(np.log(sample_scales) * 2.0),
                "type": "line",
                "name": "smoothly-varying <br>B-factor correction",
                "xaxis": "x",
                "yaxis": "y",
            }
        )
        data.append(
            {
                "x": list(smoother_phis),
                "y": list(parameters),
                "type": "scatter",
                "mode": "markers",
                "name": "smoothly-varying <br>B-factor parameters",
                "xaxis": "x",
                "yaxis": "y",
            }
        )
        if parameter_esds:
            data[-1]["error_y"] = {"type": "data", "array": list(parameter_esds)}
    d["smooth_scale_model"]["data"].extend(data)
    return d


absorption_help_msg = """
The absorption correction uses a set of spherical harmonic functions as the
basis of a smoothly varying absorption correction as a function of phi and
theta (relative to the crystal reference frame). The correction is given by:
S(s0, s1) = 1 + sum [Clm * (Ylm(s1) + Ylm(s0)) /2]
where each Ylm is a spherical harmonic and each Clm is a model parameter.
The parameters (Clm values) are shown in the 'Absorption correction surface
parameters' plot. For each l=1,2,...,lmax, there are 2l+1 m-values
from -l, -l+1,...,l+1,l, shown in the plot from left to right for each l.
s0 and s1 are the incoming and scattered beam vectors.
"""


def plot_absorption_parameters(physical_model):
    """Make a simple plot showing the absorption parameters and errors."""
    params = physical_model.components["absorption"].parameters
    param_esds = physical_model.components["absorption"].parameter_esds
    d = {
        "absorption_parameters": {
            "data": [
                {
                    "x": [i + 0.5 for i in range(len(params))],
                    "y": list(params),
                    "type": "scatter",
                    "name": "absorption parameters",
                    "xaxis": "x",
                    "yaxis": "y",
                    "mode": "markers",
                }
            ],
            "layout": {
                "title": "Absorption correction surface parameters",
                "xaxis": {"domain": [0, 1], "anchor": "y", "title": "", "tickvals": []},
                "yaxis": {"domain": [0, 1], "anchor": "x", "title": "Parameter value"},
            },
            "help": absorption_help_msg,
        }
    }
    if param_esds:
        d["absorption_parameters"]["data"][-1]["error_y"] = {
            "type": "data",
            "array": list(param_esds),
        }

    light_grey = "#d3d3d3"
    grey = "#808080"
    shapes = []
    lmax = int(-1 + (1 + len(params)) ** 0.5)
    ls = [i + 1 for i in range(lmax)]
    ns = [(2 * l) + 1 for l in ls]
    annotations = []
    start = 0
    for i, n in enumerate(ns):
        fillcolor = [light_grey, grey][i % 2]  # alternate colours
        shapes.append(
            {
                "type": "rect",
                "xref": "x",
                "yref": "paper",
                "x0": start,
                "y0": 0,
                "x1": start + n,
                "y1": 1,
                "fillcolor": fillcolor,
                "opacity": 0.2,
                "line": {"width": 0},
            }
        )
        annotations.append(
            {
                "xref": "x",
                "yref": "paper",
                "x": start + (n / 2.0),
                "y": 1,
                "text": "l=%s" % ls[i],
                "showarrow": False,
                "yshift": 20,
            }
        )
        start += n
    d["absorption_parameters"]["layout"]["shapes"] = shapes
    d["absorption_parameters"]["layout"]["annotations"] = annotations
    return d


def plot_absorption_surface(physical_model):
    """Plot an absorption surface for a physical scaling model."""

    d = {
        "absorption_surface": {
            "data": [],
            "layout": {
                "title": "Absorption correction surface",
                "xaxis": {"domain": [0, 1], "anchor": "y", "title": "theta (degrees)"},
                "yaxis": {"domain": [0, 1], "anchor": "x", "title": "phi (degrees)"},
            },
            "help": absorption_help_msg,
        }
    }

    params = physical_model.components["absorption"].parameters

    order = int(-1.0 + ((1.0 + len(params)) ** 0.5))
    lfg = scitbxmath.log_factorial_generator(2 * order + 1)
    STEPS = 50
    phi = np.linspace(0, 2 * np.pi, 2 * STEPS)
    theta = np.linspace(0, np.pi, STEPS)
    THETA, _ = np.meshgrid(theta, phi)
    lmax = int(-1.0 + ((1.0 + len(params)) ** 0.5))
    Intensity = np.ones(THETA.shape)
    counter = 0
    sqrt2 = pymath.sqrt(2)
    nsssphe = scitbxmath.nss_spherical_harmonics(order, 50000, lfg)
    for l in range(1, lmax + 1):
        for m in range(-l, l + 1):
            for it, t in enumerate(theta):
                for ip, p in enumerate(phi):
                    Ylm = nsssphe.spherical_harmonic(l, abs(m), t, p)
                    if m < 0:
                        r = sqrt2 * ((-1) ** m) * Ylm.imag
                    elif m == 0:
                        assert Ylm.imag == 0.0
                        r = Ylm.real
                    else:
                        r = sqrt2 * ((-1) ** m) * Ylm.real
                    Intensity[ip, it] += params[counter] * r
            counter += 1
    d["absorption_surface"]["data"].append(
        {
            "x": list(theta * 180.0 / np.pi),
            "y": list(phi * 180.0 / np.pi),
            "z": list(Intensity.T.tolist()),
            "type": "heatmap",
            "colorscale": "Viridis",
            "colorbar": {"title": "inverse <br>scale factor"},
            "name": "absorption surface",
            "xaxis": "x",
            "yaxis": "y",
        }
    )
    return d


def plot_outliers(data):
    """plots positions of outliers"""

    if not data["z"]:
        return {"outlier_xy_positions": {}, "outliers_vs_z": {}}

    hist = flex.histogram(
        flex.double(data["z"]), n_slots=min(100, int(len(data["z"]) * 10))
    )

    d = {
        "outlier_xy_positions": {
            "data": [
                {
                    "x": data["x"],
                    "y": data["y"],
                    "type": "scatter",
                    "mode": "markers",
                    "xaxis": "x",
                    "yaxis": "y",
                }
            ],
            "layout": {
                "title": "Outlier x-y positions",
                "xaxis": {
                    "anchor": "y",
                    "title": "x (px)",
                    "range": [0, data["image_size"][0]],
                },
                "yaxis": {
                    "anchor": "x",
                    "title": "y (px)",
                    "range": [0, data["image_size"][1]],
                },
            },
        },
        "outliers_vs_z": {
            "data": [
                {
                    "x": list(hist.slot_centers()),
                    "y": list(hist.slots()),
                    "type": "bar",
                    "name": "outliers vs rotation",
                }
            ],
            "layout": {
                "title": "Outlier distribution across frames",
                "xaxis": {"title": "frame"},
                "yaxis": {"title": "count"},
                "bargap": 0,
            },
        },
    }

    return d


def normal_probability_plot(data):
    """Plot the distribution of normal probabilities of errors."""
    norm = distributions.normal_distribution()

    n = len(data["delta_hl"])
    if n <= 10:
        a = 3 / 8
    else:
        a = 0.5

    y = flex.sorted(flex.double(data["delta_hl"]))
    x = [norm.quantile((i + 1 - a) / (n + 1 - (2 * a))) for i in range(n)]

    H, xedges, yedges = np.histogram2d(np.array(x), y.as_numpy_array(), bins=(200, 200))
    nonzeros = np.nonzero(H)
    z = np.empty(H.shape)
    z[:] = np.NAN
    z[nonzeros] = H[nonzeros]

    # also make a histogram
    histy = flex.histogram(y, n_slots=100)
    # make a gaussian for reference also
    n = y.size()
    width = histy.slot_centers()[1] - histy.slot_centers()[0]
    gaussian = []
    from math import exp, pi

    for x in histy.slot_centers():
        gaussian.append(n * width * exp(-(x ** 2) / 2.0) / ((2.0 * pi) ** 0.5))

    return {
        "normal_distribution_plot": {
            "data": [
                {
                    "x": xedges.tolist(),
                    "y": yedges.tolist(),
                    "z": z.transpose().tolist(),
                    "type": "heatmap",
                    "name": "normalised deviations",
                    "colorbar": {
                        "title": "Number of reflections",
                        "titleside": "right",
                    },
                    "colorscale": "Jet",
                },
                {
                    "x": [-5, 5],
                    "y": [-5, 5],
                    "type": "scatter",
                    "mode": "lines",
                    "name": "z = m",
                    "color": "rgb(0,0,0)",
                },
            ],
            "layout": {
                "title": "Normal probability plot with error model applied",
                "xaxis": {"anchor": "y", "title": "Order statistic medians, m"},
                "yaxis": {"anchor": "x", "title": "Ordered responses, z"},
            },
            "help": """\
This plot shows the normalised devations (of each reflection from the
group-weighted mean), sorted in order and plotted against the expected order
based on a normal distribution model. A true normal distribution of deviations
would give the straight line indicated. If the errors are well described by
this model, the ordered responses should closely fit the straight line to
high absolute values of x (>3), where there is typically a deviation away from
the line due to wide tails of the distribution.
""",
        },
        "nor_dev_hist": {
            "data": [
                {
                    "x": list(histy.slot_centers()),
                    "y": list(histy.slots()),
                    "type": "bar",
                    "name": "dataset normalised deviations",
                },
                {
                    "x": list(histy.slot_centers()),
                    "y": gaussian,
                    "type": "scatter",
                    "name": "Ideal normal distribution",
                },
            ],
            "layout": {
                "title": "Normal deviations with error model applied",
                "xaxis": {"anchor": "y", "title": "Normalised deviation"},
                "yaxis": {"anchor": "x", "title": "Number of reflections"},
            },
            "help": """\
This plot shows the distribution of normalised devations (of each reflection
from the group-weighted mean), for the reflections used to minimise the error
model. A true normal distribution is indicated.
""",
        },
    }
