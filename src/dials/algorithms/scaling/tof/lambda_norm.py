"""Wavelength normalisation of Laue/ToF intensities by Gaussian-process regression."""

from __future__ import annotations

import logging
from dataclasses import dataclass, field

import numpy as np

from cctbx import crystal, miller
from dxtbx import flumpy

from dials.algorithms.scaling.scaling_utilities import calculate_prescaling_correction
from dials.array_family import flex
from dials.util import tabulate
from dials_tof_scaling_ext import (
    GPRWavelengthScale,
    LambdaNormVarianceModel,
    leave_one_out_samples,
)

logger = logging.getLogger("dials.algorithms.scaling.tof.lambda_norm")

WAVELENGTH = "wavelength_cal"
SCALE_COLUMN = "lambda_norm.scale"
SIGMA_COLUMN = "lambda_norm.sigma"


@dataclass
class Observations:
    """Observations used in the fit, sorted by symmetry-equivalent group."""

    intensity: np.ndarray
    variance: np.ndarray
    wavelength: np.ndarray
    experiment: np.ndarray
    group_offsets: np.ndarray


@dataclass
class LambdaNormResult:
    gpr: GPRWavelengthScale
    variance_model: LambdaNormVarianceModel
    experiment_factors: np.ndarray
    cycles: list[dict] = field(default_factory=list)


def intensity_kind(reflections, choice):
    if choice == "auto":
        return "prf" if "intensity.prf.value" in reflections else "sum"
    return {"profile": "prf", "sum": "sum"}[choice]


def select_observations(reflections, experiments, params) -> Observations:
    """The observations the fit is trained on, corrected as dials.scale corrects them."""
    kind = intensity_kind(reflections, params.intensity_choice)
    flags = reflections.flags
    sel = reflections.get_flags(
        flags.integrated_prf if kind == "prf" else flags.integrated_sum
    )
    for excluded in (
        flags.user_excluded_in_scaling,
        flags.excluded_for_scaling,
        flags.bad_for_scaling,
    ):
        sel &= ~reflections.get_flags(excluded)
    sel &= reflections["miller_index"].as_vec3_double().norms() > 0
    sel &= reflections[f"intensity.{kind}.variance"] > 0
    if "partiality" in reflections:
        sel &= reflections["partiality"] >= params.partiality_cutoff
        sel &= reflections["partiality"] > 0
    if params.wavelength_range:
        lam = reflections[WAVELENGTH]
        sel &= (lam >= params.wavelength_range[0]) & (lam <= params.wavelength_range[1])
    table = calculate_prescaling_correction(reflections.select(sel))

    conversion = table["prescaling_correction"]
    if kind == "sum" and "partiality" in table:
        conversion /= table["partiality"]
    intensity = table[f"intensity.{kind}.value"] * conversion
    variance = table[f"intensity.{kind}.variance"] * conversion * conversion
    keep = intensity / flex.sqrt(variance) >= params.min_isigi
    table = table.select(keep)
    intensity = flumpy.to_numpy(intensity.select(keep))
    variance = flumpy.to_numpy(variance.select(keep))

    order, offsets = _group_by_asu(table, experiments, params.anomalous)
    return Observations(
        intensity=intensity[order],
        variance=variance[order],
        wavelength=flumpy.to_numpy(table[WAVELENGTH])[order],
        experiment=flumpy.to_numpy(table["id"])[order],
        group_offsets=offsets,
    )


def _group_by_asu(table, experiments, anomalous):
    space_groups = {e.crystal.get_space_group().type().number() for e in experiments}
    if len(space_groups) != 1:
        raise ValueError(
            "All experiments must share one space group to group equivalents"
        )
    first = experiments[0].crystal
    symmetry = crystal.symmetry(
        unit_cell=first.get_unit_cell(), space_group=first.get_space_group()
    )
    asu = (
        miller.set(symmetry, table["miller_index"], anomalous_flag=anomalous)
        .map_to_asu()
        .indices()
    )
    hkl = flumpy.to_numpy(asu.as_vec3_double()).astype(np.int64)
    _, group = np.unique(hkl, axis=0, return_inverse=True)
    group = group.ravel()
    order = np.argsort(group, kind="stable")
    offsets = np.concatenate(([0], np.cumsum(np.bincount(group))))
    return order, offsets


def fit_lambda_norm(obs: Observations, params) -> LambdaNormResult:
    """Fit w(lambda), re-forming the leave-one-out samples each cycle."""
    I = flex.double(obs.intensity)
    var = flex.double(obs.variance)
    lam = flex.double(obs.wavelength)
    offsets = flex.size_t(obs.group_offsets.astype(np.uint64))

    vmodel = LambdaNormVarianceModel()
    if params.weighting == "expected":
        vmodel.fit(lam, I, var)
        _log_variance_model(vmodel)

    lambda_ref = params.reference_wavelength or float(np.median(obs.wavelength))
    lam_min, lam_max = params.wavelength_range or (0.0, 0.0)
    gpr = GPRWavelengthScale(
        nbins=params.gpr.nbins,
        lengthscale=params.gpr.lengthscale or 0.0,
        matern32=params.gpr.kernel == "matern32",
        lambda_min=lam_min,
        lambda_max=lam_max,
    )
    result = LambdaNormResult(
        gpr=gpr,
        variance_model=vmodel,
        experiment_factors=np.ones(obs.experiment.max() + 1),
    )
    w = np.ones(len(obs.intensity))
    for cycle in range(params.gpr.n_cycles):
        k = result.experiment_factors
        samples = leave_one_out_samples(
            I,
            var,
            lam,
            flex.double(w * k[obs.experiment]),
            offsets,
            vmodel,
            params.weighting == "expected",
        )
        idx = flumpy.to_numpy(samples.index).astype(np.int64)
        rho = flumpy.to_numpy(samples.rho)
        u = flumpy.to_numpy(samples.u)
        if params.per_experiment_factor:
            k = _experiment_factors(obs.experiment[idx], rho / w[idx], u, len(k))
            result.experiment_factors = k
            e = obs.experiment[idx]
            rho = rho / k[e]
            u = u * k[e] ** 2
        ok = gpr.fit(
            flex.double(obs.wavelength[idx]),
            flex.double(rho),
            flex.double(u),
            lambda_ref,
        )
        if not ok:
            raise ValueError(f"Wavelength normalisation failed: {gpr.message()}")
        w = flumpy.to_numpy(gpr.scale(lam))
        result.cycles.append(_cycle_summary(gpr, len(idx)))
        logger.info(
            "Cycle %d: %d samples, length scale %.4f (%.3f Å), sigma_f %.4f, "
            "noise inflation %.2f, log marginal likelihood %.2f",
            cycle + 1,
            len(idx),
            gpr.lengthscale(),
            gpr.lengthscale() * gpr.lambda_ref(),
            gpr.sigma_f(),
            gpr.noise_inflation(),
            gpr.log_marginal_likelihood(),
        )
    _log_warnings(gpr)
    return result


def _experiment_factors(experiment, ratio, u, n):
    num = np.bincount(experiment, weights=u * ratio, minlength=n)
    den = np.bincount(experiment, weights=u, minlength=n)
    k = np.ones(n)
    present = den > 0
    k[present] = num[present] / den[present]
    k[~(k > 0)] = 1.0
    # w and k are degenerate in overall scale; w carries it via lambda_ref
    return k / np.exp(np.mean(np.log(k[present])))


def _cycle_summary(gpr, n_samples):
    return {
        "n_samples": n_samples,
        "lengthscale": gpr.lengthscale(),
        "sigma_f": gpr.sigma_f(),
        "noise_inflation": gpr.noise_inflation(),
        "log_marginal_likelihood": gpr.log_marginal_likelihood(),
    }


def _log_variance_model(vmodel):
    edges = list(vmodel.lambda_edges())
    if not edges:
        logger.info("Too few observations for the variance model; using measured σ")
        return
    rows = [
        [f"{lo:.3f} - {hi:.3f}", f"{b:.4g}", f"{share:.2f}"]
        for lo, hi, b, share in zip(
            edges[:-1], edges[1:], vmodel.b_values(), vmodel.counting_share()
        )
    ]
    logger.info(
        "Expected-σ weighting, σ²_exp = (σ² - b·I) + b·I_pred:\n%s",
        tabulate(rows, ["λ range (Å)", "b", "median b·I/σ²"]),
    )


def _log_warnings(gpr):
    if gpr.lengthscale_at_lower_limit():
        logger.warning(
            "Length scale is at its lower search limit; the fit may be following "
            "noise. Consider fewer gpr.nbins."
        )
    if gpr.lengthscale_at_upper_limit():
        logger.warning(
            "Length scale is at its upper search limit; w(λ) is nearly flat."
        )
    if gpr.noise_inflation_at_upper_limit():
        logger.warning(
            "Noise inflation is at its upper search limit; bin scatter is much "
            "larger than the σ predict."
        )


def apply_lambda_norm(reflections, gpr):
    """Divide intensities by w(lambda) and record w and its relative SD."""
    lam = reflections[WAVELENGTH]
    w = gpr.scale(lam)
    for kind in ("prf", "sum"):
        if f"intensity.{kind}.value" in reflections:
            reflections[f"intensity.{kind}.value"] /= w
            reflections[f"intensity.{kind}.variance"] /= w * w
    reflections[SCALE_COLUMN] = w
    reflections[SIGMA_COLUMN] = gpr.uncertainty(lam)
    return reflections


def normalisation_table(gpr, n_points=12):
    lam = np.linspace(gpr.lambda_min(), gpr.lambda_max(), n_points)
    w = gpr.scale(flex.double(lam))
    sd = gpr.uncertainty(flex.double(lam))
    rows = [[f"{l:.3f}", f"{s:.4f}", f"{d:.4f}"] for l, s, d in zip(lam, w, sd)]
    return tabulate(rows, ["λ (Å)", "w(λ)", "σ(w)/w"])


def ascii_spectrum(gpr, width=60, height=15):
    """w(λ) as a text plot, with ':' marking the reference wavelength."""
    lam_min, lam_max, ref = gpr.lambda_min(), gpr.lambda_max(), gpr.lambda_ref()
    w = np.array(gpr.scale(flex.double(np.linspace(lam_min, lam_max, width))))
    lo, hi = w.min(), w.max()
    if hi - lo < 1e-6:
        lo, hi = lo - 0.5, hi + 0.5
    grid = [[" "] * width for _ in range(height)]
    refcol = None
    if lam_min <= ref <= lam_max:
        refcol = round((ref - lam_min) / (lam_max - lam_min) * (width - 1))
        for row in grid:
            row[refcol] = ":"
    rows = np.rint((hi - w) / (hi - lo) * (height - 1)).astype(int)
    for col, row in enumerate(rows):
        grid[row][col] = "*"

    labels = np.linspace(hi, lo, height)
    lines = [f"{label:8.4f} |{''.join(row)}" for label, row in zip(labels, grid)]
    lines.append(" " * 9 + "+" + "-" * width)
    lo_label, hi_label = f"{lam_min:.3f}", f"{lam_max:.3f}"
    axis = lo_label + hi_label.rjust(width - len(lo_label))
    lines.append(" " * 10 + axis + "  λ (Å)")
    if refcol is not None:
        lines.append(" " * 10 + f"(':' marks the reference wavelength {ref:.3f} Å)")
    return "\n".join(lines)


def lambda_norm_dict(result: LambdaNormResult, n_points=201):
    """The fit and its inputs, for json output and plotting."""
    gpr = result.gpr
    lam = np.linspace(gpr.lambda_min(), gpr.lambda_max(), n_points)
    return {
        "lambda_ref": gpr.lambda_ref(),
        "lambda_min": gpr.lambda_min(),
        "lambda_max": gpr.lambda_max(),
        "lengthscale_log_lambda": gpr.lengthscale(),
        "sigma_f": gpr.sigma_f(),
        "noise_inflation": gpr.noise_inflation(),
        "log_marginal_likelihood": gpr.log_marginal_likelihood(),
        "experiment_factors": list(result.experiment_factors),
        "cycles": result.cycles,
        "variance_model": {
            "lambda_edges": list(result.variance_model.lambda_edges()),
            "b": list(result.variance_model.b_values()),
        },
        "curve": {
            "lambda": list(lam),
            "w": list(gpr.scale(flex.double(lam))),
            "sigma": list(gpr.uncertainty(flex.double(lam))),
        },
        "training": {
            "lambda": list(gpr.training_lambda()),
            "w": list(gpr.training_scale()),
            "sigma": list(gpr.training_sigma()),
        },
    }


def lambda_norm_plots(data):
    curve, bins = data["curve"], data["training"]
    w = np.array(curve["w"])
    sd = np.array(curve["sigma"])
    bw = np.array(bins["w"])
    return {
        "lambda_norm": {
            "data": [
                {
                    "x": curve["lambda"] + curve["lambda"][::-1],
                    "y": list(w * np.exp(sd)) + list((w * np.exp(-sd))[::-1]),
                    "type": "scatter",
                    "fill": "toself",
                    "line": {"width": 0},
                    "opacity": 0.3,
                    "hoverinfo": "skip",
                    "name": "±1σ",
                },
                {
                    "x": curve["lambda"],
                    "y": list(w),
                    "type": "scatter",
                    "mode": "lines",
                    "name": "w(λ)",
                },
                {
                    "x": bins["lambda"],
                    "y": list(bw),
                    "error_y": {"type": "data", "array": list(bw * bins["sigma"])},
                    "type": "scatter",
                    "mode": "markers",
                    "name": "training bins",
                },
            ],
            "layout": {
                "title": "Wavelength normalisation",
                "xaxis": {"title": "λ (Å)"},
                "yaxis": {"title": "w(λ)"},
            },
            "help": "The fitted w(λ), normalised to 1 at the reference wavelength, "
            "with the binned leave-one-out ratios it was fitted to. The curve "
            "should follow the bins without chasing individual points.",
        }
    }
