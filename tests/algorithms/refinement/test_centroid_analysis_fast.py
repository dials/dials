"""Tests for the numpy-accelerated CentroidAnalyser implementation.

Verifies that the numpy.histogram-based block counting and numpy.digitize +
numpy.bincount-based residual averaging produce results identical to the
reference per-block boolean selection implementation.
"""

from __future__ import annotations

import math
import time

import numpy as np

from dials.algorithms.refinement.analysis.centroid_analysis import CentroidAnalyser
from dials.array_family import flex

RAD2DEG = 180.0 / math.pi
DEG2RAD = math.pi / 180.0


def _reference_nr(phi_obs_deg, phi_range, nblocks, block_size):
    """Reference implementation: per-block boolean selection for counting."""
    nr = flex.int()
    for i in range(nblocks - 1):
        blk_start = phi_range[0] + i * block_size
        blk_end = blk_start + block_size
        sel = (phi_obs_deg >= blk_start) & (phi_obs_deg < blk_end)
        nr.append(sel.count(True))
    blk_start = phi_range[0] + (nblocks - 1) * block_size
    blk_end = phi_range[1]
    sel = (phi_obs_deg >= blk_start) & (phi_obs_deg <= blk_end)
    nr.append(sel.count(True))
    return nr


def _reference_block_means(phi_obs_deg, resid, phi_range, nblocks, block_size):
    """Reference implementation: per-block boolean selection for mean residuals."""
    result = flex.double()
    for i in range(nblocks - 1):
        blk_start = phi_range[0] + i * block_size
        blk_end = blk_start + block_size
        sel = (phi_obs_deg >= blk_start) & (phi_obs_deg < blk_end)
        result.append(flex.mean(resid.select(sel)))
    blk_start = phi_range[0] + (nblocks - 1) * block_size
    blk_end = phi_range[1]
    sel = (phi_obs_deg >= blk_start) & (phi_obs_deg <= blk_end)
    result.append(flex.mean(resid.select(sel)))
    return result


def _make_reflections(n, phi_min_deg=0.0, phi_max_deg=90.0, seed=42, exp_id=0):
    """Build a minimal synthetic reflections table for CentroidAnalyser."""
    rng = np.random.default_rng(seed)

    phi_deg = rng.uniform(phi_min_deg, phi_max_deg, n)
    phi_rad = phi_deg * DEG2RAD

    # Small random x/y positions (not 0,0 so they survive the (x==0)&(y==0) filter)
    x = rng.uniform(0.1, 10.0, n)
    y = rng.uniform(0.1, 10.0, n)

    # Synthetic residuals (small random offsets)
    x_resid = rng.normal(0.0, 0.01, n)
    y_resid = rng.normal(0.0, 0.01, n)
    phi_resid_rad = rng.normal(0.0, 0.001, n)

    refl = flex.reflection_table()
    refl["id"] = flex.int(np.full(n, exp_id, dtype=np.int32).tolist())
    refl["miller_index"] = flex.miller_index(
        [(1, 0, 0)] * n  # non-zero miller index to survive filter
    )
    # xyzcal.mm: x, y, phi_cal (we store phi_obs + residual as phi_cal)
    phi_cal_rad = phi_rad + phi_resid_rad
    refl["xyzcal.mm"] = flex.vec3_double(
        list(zip(x.tolist(), y.tolist(), phi_cal_rad.tolist()))
    )
    # xyzobs.mm.value: x_obs, y_obs, phi_obs
    x_obs = x - x_resid
    y_obs = y - y_resid
    refl["xyzobs.mm.value"] = flex.vec3_double(
        list(zip(x_obs.tolist(), y_obs.tolist(), phi_rad.tolist()))
    )
    # Pre-supply residuals so CentroidAnalyser skips its own computation
    refl["x_resid"] = flex.double(x_resid.tolist())
    refl["y_resid"] = flex.double(y_resid.tolist())
    refl["phi_resid"] = flex.double(phi_resid_rad.tolist())

    return refl


# ---------------------------------------------------------------------------
# Tests for __init__ block counting (numpy.histogram path)
# ---------------------------------------------------------------------------


def test_block_counting_matches_reference():
    """nref_per_block from CentroidAnalyser must match the reference implementation."""
    n = 5000
    refl = _make_reflections(n)
    ca = CentroidAnalyser(refl)
    results = ca(calc_average_residuals=False, calc_periodograms=False)

    assert len(results) == 1
    res = results[0]
    nblocks = res["nblocks"]
    block_size = res["block_size"]
    phi_range = res["phi_range"]

    # Reconstruct phi_obs_deg to run reference
    ref_this_exp = ca._reflections.select(ca._reflections["id"] == 0)
    phi_obs_deg = ref_this_exp["xyzobs.mm.value"].parts()[2] * RAD2DEG

    ref_nr = _reference_nr(phi_obs_deg, phi_range, nblocks, block_size)
    got_nr = res["nref_per_block"]

    assert list(got_nr) == list(ref_nr), (
        f"Block counts differ:\n  got={list(got_nr)}\n  ref={list(ref_nr)}"
    )


def test_block_counting_total_equals_nref():
    """Total reflections across blocks must equal the number in the experiment."""
    n = 5000
    refl = _make_reflections(n)
    ca = CentroidAnalyser(refl)
    results = ca(calc_average_residuals=False, calc_periodograms=False)
    res = results[0]
    assert flex.sum(res["nref_per_block"]) == len(ca._reflections)


# ---------------------------------------------------------------------------
# Tests for __call__ residual averaging (numpy.digitize + bincount path)
# ---------------------------------------------------------------------------


def test_residual_averages_match_reference():
    """Per-block residual means from the fast path must match the reference."""
    n = 5000
    refl = _make_reflections(n)
    ca = CentroidAnalyser(refl)
    results = ca(calc_average_residuals=True, calc_periodograms=False)

    res = results[0]
    nblocks = res["nblocks"]
    block_size = res["block_size"]
    phi_range = res["phi_range"]

    ref_this_exp = ca._reflections.select(ca._reflections["id"] == 0)
    phi_obs_deg = ref_this_exp["xyzobs.mm.value"].parts()[2] * RAD2DEG
    x_resid = ref_this_exp["x_resid"]
    y_resid = ref_this_exp["y_resid"]
    phi_resid = ref_this_exp["phi_resid"]

    ref_xr = _reference_block_means(
        phi_obs_deg, x_resid, phi_range, nblocks, block_size
    )
    ref_yr = _reference_block_means(
        phi_obs_deg, y_resid, phi_range, nblocks, block_size
    )
    ref_pr = _reference_block_means(
        phi_obs_deg, phi_resid, phi_range, nblocks, block_size
    )

    # Apply the same edge correction the implementation applies
    if nblocks > 2 and block_size < 3.0:
        ref_xr[0] = ref_xr[1]
        ref_xr[-1] = ref_xr[-2]
        ref_yr[0] = ref_yr[1]
        ref_yr[-1] = ref_yr[-2]
        ref_pr[0] = ref_pr[1]
        ref_pr[-1] = ref_pr[-2]

    got_xr = res["av_x_resid_per_block"]
    got_yr = res["av_y_resid_per_block"]
    got_pr = res["av_phi_resid_per_block"]

    tol = 1e-10
    for i, (g, r) in enumerate(zip(got_xr, ref_xr)):
        assert abs(g - r) < tol, f"av_x_resid_per_block[{i}]: got {g}, ref {r}"
    for i, (g, r) in enumerate(zip(got_yr, ref_yr)):
        assert abs(g - r) < tol, f"av_y_resid_per_block[{i}]: got {g}, ref {r}"
    for i, (g, r) in enumerate(zip(got_pr, ref_pr)):
        assert abs(g - r) < tol, f"av_phi_resid_per_block[{i}]: got {g}, ref {r}"


# ---------------------------------------------------------------------------
# Edge case: nblocks=1 (very sparse data forces single block)
# ---------------------------------------------------------------------------


def test_single_block_edge_case():
    """With very few reflections, CentroidAnalyser converges to nblocks=1."""
    # 3 reflections spread over 90 degrees — guaranteed to end up as nblocks=1
    # because min_nr < 50 and block halving will hit nblocks=1
    refl = _make_reflections(3, phi_min_deg=0.0, phi_max_deg=90.0, seed=7)
    ca = CentroidAnalyser(refl)
    results = ca(calc_average_residuals=False, calc_periodograms=False)
    assert results[0]["nblocks"] == 1
    assert flex.sum(results[0]["nref_per_block"]) == len(ca._reflections)


# ---------------------------------------------------------------------------
# Edge case: nblocks=2
# ---------------------------------------------------------------------------


def test_two_block_edge_case():
    """nblocks=2 exercised with exactly 100 reflections in each half."""
    # Place 100 reflections in [0, 45) and 100 in [45, 90] to ensure nblocks >= 2
    # with min_nr >= 50
    n_half = 100
    rng = np.random.default_rng(99)
    phi_lo = rng.uniform(0.0, 44.99, n_half)
    phi_hi = rng.uniform(45.01, 90.0, n_half)
    phi_deg = np.concatenate([phi_lo, phi_hi])
    phi_rad = phi_deg * DEG2RAD
    n = 2 * n_half

    x = rng.uniform(0.1, 10.0, n)
    y = rng.uniform(0.1, 10.0, n)
    x_resid = rng.normal(0.0, 0.01, n)
    y_resid = rng.normal(0.0, 0.01, n)
    phi_resid_rad = rng.normal(0.0, 0.001, n)

    refl = flex.reflection_table()
    refl["id"] = flex.int([0] * n)
    refl["miller_index"] = flex.miller_index([(1, 0, 0)] * n)
    phi_cal_rad = phi_rad + phi_resid_rad
    refl["xyzcal.mm"] = flex.vec3_double(
        list(zip(x.tolist(), y.tolist(), phi_cal_rad.tolist()))
    )
    x_obs = x - x_resid
    y_obs = y - y_resid
    refl["xyzobs.mm.value"] = flex.vec3_double(
        list(zip(x_obs.tolist(), y_obs.tolist(), phi_rad.tolist()))
    )
    refl["x_resid"] = flex.double(x_resid.tolist())
    refl["y_resid"] = flex.double(y_resid.tolist())
    refl["phi_resid"] = flex.double(phi_resid_rad.tolist())

    ca = CentroidAnalyser(refl)
    results = ca(calc_average_residuals=True, calc_periodograms=False)
    res = results[0]
    assert res["nblocks"] >= 2
    assert flex.sum(res["nref_per_block"]) == len(ca._reflections)
    assert len(res["av_x_resid_per_block"]) == res["nblocks"]


# ---------------------------------------------------------------------------
# Custom av_callback falls back to original loop
# ---------------------------------------------------------------------------


def test_custom_av_callback_fallback():
    """A non-flex.mean callback must still produce valid per-block results."""

    def my_median(arr):
        return flex.median(arr)

    n = 500
    refl = _make_reflections(n)
    ca = CentroidAnalyser(refl, av_callback=my_median)
    results = ca(calc_average_residuals=True, calc_periodograms=False)
    res = results[0]
    assert len(res["av_x_resid_per_block"]) == res["nblocks"]
    assert len(res["av_y_resid_per_block"]) == res["nblocks"]
    assert len(res["av_phi_resid_per_block"]) == res["nblocks"]


# ---------------------------------------------------------------------------
# Performance: __init__ must complete quickly for large datasets
# ---------------------------------------------------------------------------


def test_init_performance():
    """CentroidAnalyser.__init__ must run in under 0.3s per call for n=960000."""
    n = 960_000
    refl = _make_reflections(n, seed=1234)

    t0 = time.perf_counter()
    ca = CentroidAnalyser(refl)
    elapsed = time.perf_counter() - t0

    assert elapsed < 0.3, (
        f"CentroidAnalyser.__init__ took {elapsed:.3f}s, expected < 0.3s"
    )
    # Sanity check
    results = ca(calc_average_residuals=False, calc_periodograms=False)
    assert flex.sum(results[0]["nref_per_block"]) == len(ca._reflections)
