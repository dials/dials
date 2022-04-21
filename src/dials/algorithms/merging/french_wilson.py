from __future__ import annotations

import logging
import math
from dataclasses import dataclass

import numpy as np
from scipy import interpolate, stats

from cctbx import miller
from cctbx.french_wilson import (
    ac_zf,
    ac_zf_sd,
    ac_zj,
    ac_zj_sd,
    c_zf,
    c_zf_sd,
    c_zj,
    c_zj_sd,
)
from dxtbx import flumpy

logger = logging.getLogger(__name__)


def french_wilson(
    merged_intensities: miller.array,
    n_bins: int = 50,
    sigma_iobs_rejection_criterion: float = -4.0,
) -> miller.array:
    assert merged_intensities.is_xray_intensity_array()
    assert merged_intensities.sigmas() is not None
    assert merged_intensities.is_unique_set_under_symmetry()

    intensities = flumpy.to_numpy(merged_intensities.data())
    sigmas = flumpy.to_numpy(merged_intensities.sigmas())
    is_centric = flumpy.to_numpy(merged_intensities.centric_flags().data())
    is_acentric = ~is_centric

    rejected = (sigmas <= 0) & ((intensities <= 0) | (sigmas < 0))

    # Compute expected intensities
    four_stol_sq = 4 * flumpy.to_numpy(
        merged_intensities.sin_theta_over_lambda_sq().data()
    )
    expected_intensities = compute_expected_intensities(
        four_stol_sq,
        intensities,
        n_bins=n_bins,
        m_estimator=standardized_median,
    )
    rejected |= expected_intensities == 0
    valid = ~rejected

    J = np.zeros_like(intensities)
    sigJ = np.zeros_like(intensities)
    F = np.zeros_like(intensities)
    sigF = np.zeros_like(intensities)

    # Compute acentric moments
    posterior_moments_acentric = compute_posterior_moments_acentric(
        intensities[is_acentric], sigmas[is_acentric], expected_intensities[is_acentric]
    )
    J[is_acentric] = posterior_moments_acentric.J
    sigJ[is_acentric] = posterior_moments_acentric.sigJ
    F[is_acentric] = posterior_moments_acentric.F
    sigF[is_acentric] = posterior_moments_acentric.sigF
    valid[is_acentric] &= posterior_moments_acentric.valid

    # Compute centric moments
    posterior_moments_centric = compute_posterior_moments_centric(
        intensities[is_centric], sigmas[is_centric], expected_intensities[is_centric]
    )
    J[is_centric] = posterior_moments_centric.J
    sigJ[is_centric] = posterior_moments_centric.sigJ
    F[is_centric] = posterior_moments_centric.F
    sigF[is_centric] = posterior_moments_centric.sigF
    valid[is_centric] &= posterior_moments_centric.valid

    logger.debug(f"Rejected {np.count_nonzero(~valid)} reflections")
    return (
        merged_intensities.customized_copy(
            data=flumpy.from_numpy(F),
            sigmas=flumpy.from_numpy(sigF),
            anomalous_flag=merged_intensities.anomalous_flag(),
        )
        .set_observation_type_xray_amplitude()
        .select(flumpy.from_numpy(valid))
    )


def standardized_median(a: np.ndarray) -> float:
    return np.median(a) / math.log(2)


def compute_expected_intensities(
    x: np.ndarray,
    intensities: np.ndarray,
    n_bins: int,
    m_estimator="mean",
) -> np.ndarray:
    # Compute the mean intensities binned according to x
    bin_means, bin_edges, binnumber = stats.binned_statistic(
        x, intensities, statistic=m_estimator, bins=n_bins
    )
    bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2

    # Interpolate expected intensities
    # f = interpolate.interp1d(bin_centers, bin_means, fill_value="extrapolate")
    f = interpolate.interp1d(
        bin_centers,
        bin_means,
        bounds_error=False,
        fill_value=(bin_means[0], bin_means[-1]),
    )
    # f = interpolate.CubicSpline(bin_centers, bin_means, bc_type="clamped")
    expected_intensities = f(x)
    return expected_intensities


@dataclass
class PosteriorMoments:
    J: np.ndarray
    sigJ: np.ndarray
    F: np.ndarray
    sigF: np.ndarray
    valid: np.ndarray


def compute_posterior_moments_acentric(
    intensities: np.ndarray,
    sigmas: np.ndarray,
    expected_intensities: np.ndarray,
    h_min: int = -4,
) -> PosteriorMoments:

    h = (intensities / sigmas) - np.abs(sigmas / expected_intensities)
    logger.debug(f"h range: {h.min()} - {h.max}")
    i_sig_min = h_min + 0.3
    valid = (intensities / sigmas >= i_sig_min) & (h >= h_min)

    J = np.empty_like(h)
    sigJ = np.empty_like(h)
    F = np.empty_like(h)
    sigF = np.empty_like(h)

    # Set up interpolation functions as a function of h
    n_interp = len(ac_zj)
    x = (np.array(range(n_interp))) / 10 - 4
    J_h = interpolate.interp1d(x, ac_zj)
    sigJ_h = interpolate.interp1d(x, ac_zj_sd)
    F_h = interpolate.interp1d(x, ac_zf)
    sigF_h = interpolate.interp1d(x, ac_zf_sd)

    # Interpolation only valid for -4 < h < 3
    sel = valid & (h < 3)
    J[sel] = J_h(h[sel]) * sigmas[sel]
    sigJ[sel] = sigJ_h(h[sel]) * sigmas[sel]
    F[sel] = F_h(h[sel]) * np.sqrt(sigmas)[sel]
    sigF[sel] = sigF_h(h[sel]) * np.sqrt(sigmas)[sel]

    # For h >= 3 approximate as a normal distribution
    sel = valid & (h >= 3)
    J[sel] = h[sel] * sigmas[sel]
    sigJ[sel] = sigmas[sel]
    F[sel] = np.sqrt(J[sel])
    sigF[sel] = 0.5 * sigmas[sel] / F[sel]

    return PosteriorMoments(
        J=J,
        sigJ=sigJ,
        F=F,
        sigF=sigF,
        valid=valid,
    )


def compute_posterior_moments_centric(
    intensities: np.ndarray,
    sigmas: np.ndarray,
    expected_intensities: np.ndarray,
    h_min: int = -4,
) -> PosteriorMoments:

    h = (intensities / sigmas) - np.abs(sigmas / (2 * expected_intensities))
    logger.debug(f"h range: {h.min()} - {h.max}")
    i_sig_min = h_min + 0.3
    valid = (intensities / sigmas >= i_sig_min) & (h >= h_min)

    J = np.empty_like(h)
    sigJ = np.empty_like(h)
    F = np.empty_like(h)
    sigF = np.empty_like(h)

    # Set up interpolation functions as a function of h
    n_interp = len(c_zj)
    xh = (np.array(range(n_interp))) / 10 - 4
    J_h = interpolate.interp1d(xh, c_zj)
    sigJ_h = interpolate.interp1d(xh, c_zj_sd)
    F_h = interpolate.interp1d(xh, c_zf)
    sigF_h = interpolate.interp1d(xh, c_zf_sd)

    # Obtain posterior moments via interpolation for -4 < h < 4
    sel = valid & (h < 4)
    J[sel] = J_h(h[sel]) * sigmas[sel]
    sigJ[sel] = sigJ_h(h[sel]) * sigmas[sel]
    F[sel] = F_h(h[sel]) * np.sqrt(sigmas)[sel]
    sigF[sel] = sigF_h(h[sel]) * np.sqrt(sigmas)[sel]

    # adapted from French-Wilson w/ added x^6 term in the expansion
    sel = valid & (h >= 4)
    h_2 = 1 / np.square(h[sel])
    h_4 = np.square(h_2)
    h_6 = h_2 * h_4
    posterior_F = np.sqrt(h[sel]) * (
        1 - (3 / 8) * h_2 - (87 / 128) * h_4 - (2889 / 1024) * h_6
    )
    posterior_sigF = np.sqrt(
        h[sel] * ((1 / 4) * h_2 + (15 / 32) * h_4 + (273 / 128) * h_6)
    )
    J[sel] = h[sel] * sigmas[sel] * (1 - (1 / 2) * h_2 - (3 / 4) * h_4 - 3 * h_6)
    sigJ[sel] = 2 * sigmas[sel] * posterior_F * posterior_sigF
    F[sel] = posterior_F * np.sqrt(sigmas[sel])
    sigF[sel] = posterior_sigF * np.sqrt(sigmas[sel])

    return PosteriorMoments(
        J=J,
        sigJ=sigJ,
        F=F,
        sigF=sigF,
        valid=valid,
    )