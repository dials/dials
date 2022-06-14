"""
Make plotly plots for html output by dials.scale, dials.report or xia2.report.
"""

from __future__ import annotations

import itertools
import math

import numpy as np
from scipy.stats import norm

from dxtbx import flumpy
from scitbx import math as scitbxmath

from dials.array_family import flex
from dials_scaling_ext import calc_lookup_index, calc_theta_phi


def _get_smooth_plotting_data_from_model(model, component="scale"):
    """Return a tuple of phis, parameters, parameter esds,
    sample positions for plotting and sample scale values."""
    valid_osc = model.configdict["valid_osc_range"]
    if model.id_ == "physical":
        sample_values = flex.double(
            np.linspace(
                valid_osc[0],
                valid_osc[1],
                int((valid_osc[1] - valid_osc[0]) / 0.1) + 1,
                endpoint=True,
            )
        )  # Make a grid of
        # points with 10 points per degree.
    elif model.id_ == "dose_decay":
        sample_values = flex.double(
            np.linspace(
                0.0,
                abs(valid_osc[1] - valid_osc[0]),
                int((abs(valid_osc[1] - valid_osc[0])) / 0.1) + 1,
                endpoint=True,
            )
        )
    if component == "scale":
        if "scale" in model.configdict["corrections"]:
            scale_SF = model.components["scale"]
            norm_vals = sample_values / model.configdict["scale_rot_interval"]
            scale_SF.data = {"x": norm_vals}
            scale_SF.update_reflection_data()
            s = scale_SF.calculate_scales()
            offset = valid_osc[0]
            if model.id_ == "dose_decay":
                offset = 0
            smoother_phis = [
                (i * model.configdict["scale_rot_interval"]) + offset
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
        if "decay" in model.configdict["corrections"]:
            decay_SF = model.components["decay"]
            norm_vals = sample_values / model.configdict["decay_rot_interval"]
            decay_SF.data = {"x": norm_vals, "d": flex.double(norm_vals.size(), 1.0)}
            decay_SF.update_reflection_data()
            s = decay_SF.calculate_scales()
            smoother_phis = [
                (i * model.configdict["decay_rot_interval"]) + valid_osc[0]
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


def plot_relative_Bs(relative_Bs):
    """Make scatter plot of relative Bs for dose-decay"""
    d = {
        "dose_relative_Bs": {
            "data": [
                {
                    "x": list(range(0, len(relative_Bs))),
                    "y": relative_Bs,
                    "type": "scatter",
                    "name": "Relative B-factor per dataset",
                    "mode": "markers",
                }
            ],
            "layout": {
                "title": "Relative B-factor per dataset",
                "xaxis": {"anchor": "y", "title": "Sweep id"},
                "yaxis": {"anchor": "x", "title": "Relative B (Angstrom^2)"},
            },
        }
    }
    return d


def plot_dose_decay(dose_decay_model):
    """Plot the decay and scale corrections for the dose-decay model."""
    d = {
        "dose_decay": {
            "data": [],
            "layout": {
                "title": "Dose-decay model corrections",
                "xaxis": {
                    "domain": [0, 1],
                    "anchor": "y",
                    "title": "phi (degrees)",
                },
                "yaxis": {
                    "domain": [0, 0.45],
                    "anchor": "x",
                    "title": "Decay correction <br>scale factor",
                },
                "yaxis2": {
                    "domain": [0.50, 0.95],
                    "anchor": "x",
                    "title": "Scale correction <br>scale factor",
                },
            },
        }
    }

    data = []

    if "scale" in dose_decay_model.components:
        data = _add_smooth_scales_to_data(dose_decay_model, data, yaxis="y2")

    if "decay" in dose_decay_model.components:
        data = _add_decay_model_scales_to_data(
            dose_decay_model,
            data,
            yaxis="y",
            resolution=3.0,
        )

    d["dose_decay"]["data"] = data

    return d


def _add_smooth_scales_to_data(physical_model, data, yaxis="y2"):
    (
        smoother_phis,
        parameters,
        parameter_esds,
        sample_values,
        sample_scales,
    ) = _get_smooth_plotting_data_from_model(physical_model, component="scale")

    data.append(
        {
            "x": list(sample_values),
            "y": list(sample_scales),
            "type": "line",
            "name": "smoothly-varying <br>scale correction",
            "xaxis": "x",
            "yaxis": yaxis,
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
            "yaxis": yaxis,
        }
    )
    if parameter_esds:
        data[-1]["error_y"] = {"type": "data", "array": list(parameter_esds)}
    return data


def _add_decay_model_scales_to_data(model, data, yaxis="y", resolution=3.0):
    configdict = model.configdict
    valid_osc = configdict["valid_osc_range"]

    if model.id_ == "physical":
        sample_values = flex.double(
            np.linspace(
                valid_osc[0],
                valid_osc[1],
                int((valid_osc[1] - valid_osc[0]) / 0.1) + 1,
                endpoint=True,
            )
        )  # Make a grid of
        # points with 10 points per degree.
    elif model.id_ == "dose_decay":
        sample_values = flex.double(
            np.linspace(
                0.0,
                abs(valid_osc[1] - valid_osc[0]),
                int((abs(valid_osc[1] - valid_osc[0])) / 0.1) + 1,
                endpoint=True,
            )
        )
    d = flex.double(sample_values.size(), resolution)

    if "decay" in model.components:
        decay_SF = model.components["decay"]
        decay_SF.data = {"x": sample_values, "d": d}
        decay_SF.update_reflection_data()
        s = decay_SF.calculate_scales()
    data.append(
        {
            "x": list(sample_values),
            "y": list(s),
            "type": "line",
            "name": f"Decay scale factor <br>at {resolution} Angstrom",
            "xaxis": "x",
            "yaxis": yaxis,
        }
    )
    return data


def plot_smooth_scales(physical_model):
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
                    "domain": [0.0, 0.30],
                    "anchor": "x",
                    "title": "scale factor",
                },
                "yaxis2": {
                    "domain": [0.35, 0.65],
                    "anchor": "x",
                    "title": "relative <br> B-factor",
                },
                "yaxis3": {
                    "domain": [0.70, 1.0],
                    "anchor": "x",
                    "title": "scale factor",
                },
            },
            "help": smooth_help_msg,
        }
    }

    data = []

    if "scale" in physical_model.components:
        data = _add_smooth_scales_to_data(physical_model, data, yaxis="y3")
    if "decay" in physical_model.components:
        (
            smoother_phis,
            parameters,
            parameter_esds,
            sample_values,
            sample_scales,
        ) = _get_smooth_plotting_data_from_model(physical_model, component="decay")

        data.append(
            {
                "x": list(sample_values),
                "y": list(np.log(sample_scales) * 2.0),
                "type": "line",
                "name": "smoothly-varying <br>B-factor correction",
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
                "name": "smoothly-varying <br>B-factor parameters",
                "xaxis": "x",
                "yaxis": "y2",
            }
        )
        if parameter_esds:
            data[-1]["error_y"] = {"type": "data", "array": list(parameter_esds)}

        data = _add_decay_model_scales_to_data(
            physical_model, data, yaxis="y", resolution=3.0
        )
    d["smooth_scale_model"]["data"].extend(data)
    return d


absorption_help_msg = """
This plot shows the smoothly-varying absorption surface used to correct the
data for the effects of absorption. It is important to note that this plot does
not show the correction applied; the applied correction for a given reflection
with scattering vectors s0, s1 is given by the average of the two values on this
surface where the surface intersects those scattering vectors. The plot is in the
crystal reference frame, and the pole (polar angle 0) corresponds to the laboratory
x-axis.

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
                "text": f"l={ls[i]}",
                "showarrow": False,
                "yshift": 20,
            }
        )
        start += n
    d["absorption_parameters"]["layout"]["shapes"] = shapes
    d["absorption_parameters"]["layout"]["annotations"] = annotations
    return d


def plot_absorption_plots(physical_model, reflection_table=None):
    """Make a number of plots to help with the interpretation of the
    absorption correction."""
    # First plot the absorption surface

    d = {
        "absorption_surface": {
            "data": [],
            "layout": {
                "title": "Absorption correction surface",
                "xaxis": {
                    "domain": [0, 1],
                    "anchor": "y",
                    "title": "azimuthal angle (degrees)",
                },
                "yaxis": {
                    "domain": [0, 1],
                    "anchor": "x",
                    "title": "polar angle (degrees)",
                },
            },
            "help": absorption_help_msg,
        }
    }

    params = physical_model.components["absorption"].parameters

    order = int(-1.0 + ((1.0 + len(params)) ** 0.5))
    lfg = scitbxmath.log_factorial_generator(2 * order + 1)
    STEPS = 50
    azimuth_ = np.linspace(0, 2 * np.pi, 2 * STEPS)
    polar_ = np.linspace(0, np.pi, STEPS)
    THETA, _ = np.meshgrid(azimuth_, polar_, indexing="ij")
    lmax = int(-1.0 + ((1.0 + len(params)) ** 0.5))
    Intensity = np.ones(THETA.shape)
    undiffracted_intensity = np.ones(THETA.shape)
    counter = 0
    sqrt2 = math.sqrt(2)
    nsssphe = scitbxmath.nss_spherical_harmonics(order, 50000, lfg)
    for l in range(1, lmax + 1):
        for m in range(-l, l + 1):
            for it, t in enumerate(polar_):
                for ip, p in enumerate(azimuth_):
                    Ylm = nsssphe.spherical_harmonic(l, abs(m), t, p)
                    if m < 0:
                        r = sqrt2 * ((-1) ** m) * Ylm.imag
                    elif m == 0:
                        assert Ylm.imag == 0.0
                        r = Ylm.real
                    else:
                        r = sqrt2 * ((-1) ** m) * Ylm.real
                    Intensity[ip, it] += params[counter] * r
                    # for the undiffracted intensity, we want to add the correction
                    # at each point to the parity conjugate. We can use the fact
                    # that the odd l terms are parity odd, and even are even, to
                    # just calculate the even terms as follows
                    if l % 2 == 0:
                        undiffracted_intensity[ip, it] += params[counter] * r
            counter += 1
    d["absorption_surface"]["data"].append(
        {
            "x": list(azimuth_ * 180.0 / np.pi),
            "y": list(polar_ * 180.0 / np.pi),
            "z": list(Intensity.T.tolist()),
            "type": "heatmap",
            "colorscale": "Viridis",
            "colorbar": {"title": "inverse <br>scale factor"},
            "name": "absorption surface",
            "xaxis": "x",
            "yaxis": "y",
        }
    )

    d["undiffracted_absorption_surface"] = {
        "data": [],
        "layout": {
            "title": "Undiffracted absorption correction",
            "xaxis": {
                "domain": [0, 1],
                "anchor": "y",
                "title": "azimuthal angle (degrees)",
            },
            "yaxis": {
                "domain": [0, 1],
                "anchor": "x",
                "title": "polar angle (degrees)",
            },
        },
        "help": """
This plot shows the calculated relative absorption for a paths travelling
straight through the crystal at a given direction in a crystal-fixed frame of
reference (in spherical coordinates). This gives an indication of the effective
shape of the crystal for absorbing x-rays. In this plot, the pole (polar angle 0)
corresponds to the laboratory x-axis.
""",
    }

    d["undiffracted_absorption_surface"]["data"].append(
        {
            "x": list(azimuth_ * 180.0 / np.pi),
            "y": list(polar_ * 180.0 / np.pi),
            "z": list(undiffracted_intensity.T.tolist()),
            "type": "heatmap",
            "colorscale": "Viridis",
            "colorbar": {"title": "inverse <br>scale factor"},
            "name": "Undiffracted absorption correction",
            "xaxis": "x",
            "yaxis": "y",
        }
    )

    if not reflection_table:
        return d

    # now plot the directions of the scattering vectors

    d["vector_directions"] = {
        "data": [],
        "layout": {
            "title": "Scattering vectors in crystal frame",
            "xaxis": {
                "domain": [0, 1],
                "anchor": "y",
                "title": "azimuthal angle (degrees)",
                "range": [0, 360],
            },
            "yaxis": {
                "domain": [0, 1],
                "anchor": "x",
                "title": "polar angle (degrees)",
                "range": [0, 180],
            },
            "coloraxis": {
                "showscale": False,
            },
        },
        "help": """
This plot shows the scattering vector directions in the crystal reference frame
used to determine the absorption correction. The s0 vectors are plotted in yellow,
the s1 vectors are plotted in teal. This gives an indication of which parts of
the absorption correction surface are sampled when determining the absorption
correction. In this plot, the pole (polar angle 0) corresponds to the laboratory
x-axis.""",
    }

    STEPS = 180  # do one point per degree
    azimuth_ = np.linspace(0, 2 * np.pi, 2 * STEPS)
    polar_ = np.linspace(0, np.pi, STEPS)
    THETA, _ = np.meshgrid(azimuth_, polar_, indexing="ij")
    Intensity = np.full(THETA.shape, np.NAN)

    # note, the s1_lookup, s0_lookup is only calculated for large datasets, so
    # for small datasets we need to calculate again.
    if "s1_lookup" not in physical_model.components["absorption"].data:
        s1_lookup = calc_lookup_index(
            calc_theta_phi(reflection_table["s1c"]), points_per_degree=1
        )
        idx_polar, idx_azimuth = np.divmod(np.unique(s1_lookup), 360)
        Intensity[idx_azimuth, idx_polar] = 1
    else:
        s1_lookup = np.unique(physical_model.components["absorption"].data["s1_lookup"])
        # x is phi, y is theta
        idx_polar, idx_azimuth = np.divmod(s1_lookup, 720)
        idx_polar = idx_polar // 2  # convert from two points per degree to one
        idx_azimuth = idx_azimuth // 2
        Intensity[idx_azimuth, idx_polar] = 1

    d["vector_directions"]["data"].append(
        {
            "x": list(azimuth_ * 180.0 / np.pi),
            "y": list(polar_ * 180.0 / np.pi),
            "z": list(Intensity.T.tolist()),
            "type": "heatmap",
            "colorscale": "Viridis",
            "showscale": False,
            "xaxis": "x",
            "yaxis": "y",
            "zmin": 0,
            "zmax": 2,
        }
    )

    Intensity = np.full(THETA.shape, np.NAN)

    if "s0_lookup" not in physical_model.components["absorption"].data:
        s0_lookup = calc_lookup_index(
            calc_theta_phi(reflection_table["s0c"]), points_per_degree=1
        )
        idx_polar, idx_azimuth = np.divmod(np.unique(s0_lookup), 360)
        Intensity[idx_azimuth, idx_polar] = 2
    else:
        s0_lookup = np.unique(physical_model.components["absorption"].data["s0_lookup"])
        # x is phi, y is theta
        idx_polar, idx_azimuth = np.divmod(s0_lookup, 720)
        idx_polar = idx_polar // 2  # convert from two points per degree to one
        idx_azimuth = idx_azimuth // 2
        Intensity[idx_azimuth, idx_polar] = 2

    d["vector_directions"]["data"].append(
        {
            "x": list(azimuth_ * 180.0 / np.pi),
            "y": list(polar_ * 180.0 / np.pi),
            "z": list(Intensity.T.tolist()),
            "type": "heatmap",
            "colorscale": "Viridis",
            "showscale": False,
            "xaxis": "x",
            "yaxis": "y",
            "zmin": 0,
            "zmax": 2,
        }
    )

    scales = physical_model.components["absorption"].calculate_scales()
    hist = flex.histogram(scales, n_slots=min(100, int(scales.size() * 10)))

    d["absorption_corrections"] = {
        "data": [
            {
                "x": list(hist.slot_centers()),
                "y": list(hist.slots()),
                "type": "bar",
                "name": "Applied absorption corrections",
            },
        ],
        "layout": {
            "title": "Applied absorption corrections",
            "xaxis": {"anchor": "y", "title": "Inverse scale factor"},
            "yaxis": {"anchor": "x", "title": "Number of reflections"},
        },
    }

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


def error_model_variance_plot(data, label=None):
    bin_variances = data["binning_info"]["bin_variances"]
    initial_variances = data["binning_info"]["initial_variances"]
    xs = data["binning_info"]["bin_boundaries"]
    x = list(range(1, len(xs)))
    x_labels = [
        str(round(xs[i], 1)) + " - " + str(round(xs[i + 1], 1))
        for i in range(len(xs) - 1)
    ]
    key = (
        f"error_model_variances_{label}"
        if label is not None
        else "error_model_variances"
    )
    title = "Error model variances of normalised deviations"
    title = title + f" (error model {label})" if label is not None else title
    d = {
        key: {
            "data": [
                {
                    "x": x,
                    "y": list(initial_variances)[::-1],
                    "type": "scatter",
                    "mode": "markers",
                    "xaxis": "x",
                    "yaxis": "y",
                    "name": "Uncorrected variances",
                },
                {
                    "x": x,
                    "y": list(bin_variances)[::-1],
                    "type": "scatter",
                    "mode": "markers",
                    "xaxis": "x",
                    "yaxis": "y",
                    "name": "Corrected variances",
                },
                {
                    "x": [1, 10],
                    "y": [1, 1],
                    "type": "line",
                    "name": "Ideal normal distribution",
                },
            ],
            "layout": {
                "title": title,
                "xaxis": {
                    "anchor": "y",
                    "title": "Intensity range, expected unscaled intensity (counts)",
                    "tickvals": x,
                    "ticktext": x_labels[::-1],
                },
                "yaxis": {"anchor": "x", "title": "Variance of normalised deviations"},
            },
            "help": """\
This plot shows the variance of the normalised deviations for a given intensity
range for the error model minimisation. The expectation is that after a successful
error model correction, the variance will be one across the intensity range
measured i.e. the intensities are normally distributed about the symmetry group
best estimate, with a variance of one sigma.
""",
        }
    }
    return d


def error_regression_plot(data, label=None):
    """Plot the data from the regression fit."""
    x = data["regression_x"]
    y = data["regression_y"]
    fit = (x * (data["model_a"] ** 2)) + ((data["model_a"] * data["model_b"]) ** 2)
    key = f"regression_fit_{label}" if label is not None else "regression_fit"
    title = "Error model regression plot"
    title = title + f" (error model {label})" if label is not None else title
    return {
        key: {
            "data": [
                {
                    "x": list(x),
                    "y": list(y),
                    "type": "scatter",
                    "mode": "markers",
                    "name": "expected vs observed",
                },
                {
                    "x": list(x),
                    "y": list(fit),
                    "type": "scatter",
                    "name": "best least-squares fit",
                },
            ],
            "layout": {
                "title": title,
                "xaxis": {"anchor": "y", "title": "1/(I/sigma_obs) ^ 2"},
                "yaxis": {"anchor": "x", "title": "1/(I/sigma) ^ 2 "},
            },
            "help": """\
This plot shows the data used for the regression fit for error model refinement.
The form derives from the equation sigma^2(obs) = a^2[sigma^2 + (bI)^2].
Dividing by I^2 gives 1/(I/sigma_obs) ^ 2 = [a^2 * 1/(I/sigma) ^ 2] + [a^2 b^2]
i.e. y = mx + c.
Here, sigma_obs is the observed standard deviation of a group of symmetry
equivalents i.e. sigma_obs^2 = (Sum (I - g<Ih>)^2) / N-1.
""",
        }
    }


def normal_probability_plot(data, label=None):
    """Plot the distribution of normal probabilities of errors."""

    n = data["delta_hl"].size
    y = np.sort(data["delta_hl"])
    delta = 0.5 / n
    v = np.linspace(start=delta, stop=1.0 - delta, endpoint=True, num=n)
    x = norm.ppf(v)

    H, xedges, yedges = np.histogram2d(x, y, bins=(200, 200))
    nonzeros = np.nonzero(H)
    z = np.empty(H.shape)
    z[:] = np.NAN
    z[nonzeros] = H[nonzeros]

    # also make a histogram
    histy = flex.histogram(flumpy.from_numpy(y), n_slots=100)
    # make a gaussian for reference also
    n = y.size
    width = histy.slot_centers()[1] - histy.slot_centers()[0]
    gaussian = [
        n * width * math.exp(-(sc**2) / 2.0) / ((2.0 * math.pi) ** 0.5)
        for sc in histy.slot_centers()
    ]
    key = (
        f"normal_distribution_plot_{label}"
        if label is not None
        else "normal_distribution_plot"
    )
    title = "Normal probability plot with error model applied"
    title = title + f" (error model {label})" if label is not None else title
    key_hist = f"nor_dev_hist_{label}" if label is not None else "nor_dev_hist"
    title_hist = "Normal deviations with error model applied"
    title_hist = (
        title_hist + f" (error model {label})" if label is not None else title_hist
    )
    return {
        key: {
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
                    "colorscale": "Viridis",
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
                "title": title,
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
        key_hist: {
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
                "title": title_hist,
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


def _smooth_params_to_bins(n_param):
    if n_param > 4:
        return n_param - 2
    return n_param - 1


def plot_array_modulation_plot(array_model):
    modulation_comp = array_model.components["modulation"]
    configdict = array_model.configdict
    nxbins = _smooth_params_to_bins(configdict["n_x_mod_param"])
    nybins = _smooth_params_to_bins(configdict["n_y_mod_param"])
    sample_x_values = flex.double(
        np.linspace(0, nxbins, (nxbins * 10 + 1), endpoint=True)
    )
    sample_y_values = flex.double(
        np.linspace(0, nybins, (nybins * 10) + 1, endpoint=True)
    )
    x_vals = flex.double(np.tile(sample_x_values, len(sample_y_values)))
    y_vals = flex.double(np.repeat(sample_y_values, len(sample_x_values)))
    modulation_comp.data = {"x": x_vals, "y": y_vals}
    modulation_comp.update_reflection_data()
    z = modulation_comp.calculate_scales()

    xtickvals = flex.double(range(nxbins + 1))
    ytickvals = flex.double(range(nybins + 1))
    xticks = (xtickvals * configdict["x_det_bin_width"]) + configdict["xmin"]
    yticks = (ytickvals * configdict["y_det_bin_width"]) + configdict["ymin"]
    xticktext = [str(int(i)) for i in xticks]
    yticktext = [str(int(i)) for i in yticks]

    return {
        "array_modulation_plot": {
            "data": [
                {
                    "x": list(x_vals),
                    "y": list(y_vals),
                    "z": list(z),
                    "type": "heatmap",
                    "colorscale": "Viridis",
                }
            ],
            "layout": {
                "title": "Array model modulation correction",
                "xaxis": {
                    "title": "detector x-position (px)",
                    "showgrid": False,
                    "tickvals": list(xtickvals),
                    "ticktext": xticktext,
                },
                "yaxis": {
                    "title": "detector y-position (px)",
                    "anchor": "x",
                    "showgrid": False,
                    "tickvals": list(ytickvals),
                    "ticktext": yticktext,
                },
                "width": 500,
                "height": 450,
            },
        }
    }


def plot_array_absorption_plot(array_model):
    absorption_comp = array_model.components["absorption"]
    configdict = array_model.configdict
    nxbins = _smooth_params_to_bins(configdict["n_x_param"])
    nybins = _smooth_params_to_bins(configdict["n_y_param"])
    ntimebins = _smooth_params_to_bins(configdict["n_time_param"])

    sample_x_values = flex.double(np.linspace(0, nxbins, nxbins + 1, endpoint=True))
    sample_y_values = flex.double(np.linspace(0, nybins, nybins + 1, endpoint=True))
    sample_time_values = flex.double(
        np.linspace(0, ntimebins, ntimebins + 1, endpoint=True)
    )
    # the x-y values at one timepoint
    x_vals = np.tile(sample_x_values, len(sample_y_values))
    y_vals = np.repeat(sample_y_values, len(sample_x_values))
    # now tile this for n_time_vals
    x_vals = np.tile(x_vals, len(sample_time_values))
    y_vals = np.tile(y_vals, len(sample_time_values))

    z_vals = flex.double(
        np.repeat(sample_time_values, len(sample_x_values) * len(sample_y_values))
    )

    absorption_comp.data = {
        "x": flex.double(x_vals),
        "y": flex.double(y_vals),
        "z": flex.double(z_vals),
    }
    absorption_comp.update_reflection_data()
    z = absorption_comp.calculate_scales()

    xs = np.repeat(
        np.array(
            sample_time_values * configdict["time_rot_interval"]
            + configdict["valid_osc_range"][0]
        ),
        len(sample_x_values) * len(sample_y_values),
    )
    ys = list(
        itertools.chain(
            *itertools.repeat(
                range(len(sample_x_values) * len(sample_y_values)),
                len(sample_time_values),
            )
        )
    )

    return {
        "array_absorption_plot": {
            "data": [
                {
                    "x": list(xs),
                    "y": list(ys),
                    "z": list(z),
                    "type": "heatmap",
                    "colorscale": "Viridis",
                }
            ],
            "layout": {
                "title": "Array model absorption correction",
                "xaxis": {"title": "rotation angle (degrees)", "showgrid": False},
                "yaxis": {
                    "title": "detector x-y positional index",
                    "anchor": "x",
                    "showgrid": False,
                },
                "width": 500,
                "height": 450,
            },
            "help": """\
This plot shows a representation of the 3D absorption correction unwrapped into
2D. For each value of rotation (separated by ~decay_interval), the y-axis shows
the correction at the corners of a grid of n_absorption_bins x n_absorption_bins.
e.g. if n_absorption_bins = 3, there are 16 corners of the grid that spans the
detector surface, so the positional indices are the points (0, 0), (0, 1), (0, 2),
(0, 3), (1, 0), (1, 1) etc. on this grid. Although the plot is granular, the
correction applied during scaling is smoothly interpolated between the parameters
in 3D.
""",
        }
    }


def plot_array_decay_plot(array_model):

    decay_comp = array_model.components["decay"]
    configdict = array_model.configdict

    valid_osc = configdict["valid_osc_range"]
    n_points = max(int(math.ceil(valid_osc[1] - valid_osc[0])), 50)
    sample_x_values = flex.double(
        np.linspace(valid_osc[0], valid_osc[1], n_points + 1, endpoint=True)
    )
    n_y_bins = _smooth_params_to_bins(configdict["n_res_param"])
    sample_y_values = np.linspace(0, n_y_bins, int(n_y_bins / 0.1) + 1, endpoint=True)

    norm_x_vals = (
        flex.double(np.tile(sample_x_values, len(sample_y_values)))
        / configdict["time_rot_interval"]
    )
    norm_y_vals = flex.double(np.repeat(sample_y_values, len(sample_x_values)))
    decay_comp.data = {"x": norm_y_vals, "y": norm_x_vals}
    decay_comp.update_reflection_data()
    z = decay_comp.calculate_scales()
    x = norm_x_vals * configdict["time_rot_interval"]
    y = norm_y_vals

    tickvals = flex.double(range(n_y_bins + 1))
    resmin = (tickvals * configdict["res_bin_width"]) + configdict["resmin"]
    d = 1.0 / flex.sqrt((tickvals * configdict["res_bin_width"]) + resmin)
    ticktext = [f"{i:.3f}" for i in d]

    return {
        "array_decay_plot": {
            "data": [
                {
                    "x": list(x),
                    "y": list(y),
                    "z": list(z),
                    "type": "heatmap",
                    "colorscale": "Viridis",
                }
            ],
            "layout": {
                "title": "Array model decay correction",
                "xaxis": {"title": "rotation angle (degrees)", "showgrid": False},
                "yaxis": {
                    "title": "resolution/d-value (Angstrom)",
                    "anchor": "x",
                    "showgrid": False,
                    "tickvals": list(tickvals),
                    "ticktext": ticktext,
                },
                "width": 500,
                "height": 450,
            },
        }
    }
