"""Sauter-Poon outlier rejection on a detector whose panels have different
pixel sizes: residuals must be converted to mm with each panel's own pixel
size, otherwise panels with small pixels look worse than they are."""

from __future__ import annotations

import random

from dials.algorithms.refinement.outlier_detection.sauter_poon import SauterPoon
from dials.array_family import flex


def _make_reflections(px_sizes, n_per_panel=200, n_outliers=10, seed=1):
    random.seed(seed)
    refl = flex.reflection_table()
    hkl, obs, calc, panel, is_outlier = [], [], [], [], []
    for ipanel, (px_x, px_y) in enumerate(px_sizes):
        for i in range(n_per_panel):
            x, y = random.uniform(10, 500), random.uniform(10, 500)
            # 0.02 mm rms residual on every panel, expressed in that panel's pixels
            dx, dy = random.gauss(0, 0.02) / px_x, random.gauss(0, 0.02) / px_y
            outlier = i < n_outliers
            if outlier:
                dx += random.choice((-1, 1)) * 0.3 / px_x
                dy += random.choice((-1, 1)) * 0.3 / px_y
            hkl.append((ipanel, i, 0))
            obs.append((x + dx, y + dy, 0.0))
            calc.append((x, y, 0.0))
            panel.append(ipanel)
            is_outlier.append(outlier)
    refl["miller_index"] = flex.miller_index(hkl)
    refl["xyzobs.px.value"] = flex.vec3_double(obs)
    refl["xyzcal.px"] = flex.vec3_double(calc)
    refl["panel"] = flex.size_t(panel)
    refl["id"] = flex.int(len(refl), 0)
    refl.set_flags(flex.bool(len(refl), True), refl.flags.predicted)
    return refl, flex.bool(is_outlier)


def test_sauter_poon_per_panel_pixel_sizes():
    px_sizes = [(0.2, 0.2), (0.05, 0.05)]
    refl, truth = _make_reflections(px_sizes)

    # pooled over panels with per-panel pixel sizes: the inserted outliers
    # are found on both panels and few good reflections are lost
    od = SauterPoon(px_sz=px_sizes, separate_panels=False, min_num_obs=20)
    assert od(refl)
    flagged = refl.get_flags(refl.flags.centroid_outlier)
    assert (flagged & truth).count(True) == truth.count(True)
    assert (flagged & ~truth).count(True) < 0.05 * (~truth).count(True)

    # with a single (wrong for one panel) pixel size the small-pixel panel is
    # over-rejected when pooled: this is the behaviour the fix addresses
    refl2, _ = _make_reflections(px_sizes)
    od_single = SauterPoon(px_sz=px_sizes[0], separate_panels=False, min_num_obs=20)
    assert od_single(refl2)
    flagged2 = refl2.get_flags(refl2.flags.centroid_outlier)
    small_px = refl2["panel"] == 1
    assert (flagged2 & small_px).count(True) > (flagged & small_px).count(True)

    # a single pixel size still works for an ordinary detector
    refl3, truth3 = _make_reflections([(0.1, 0.1)])
    od3 = SauterPoon(px_sz=(0.1, 0.1), separate_panels=True, min_num_obs=20)
    assert od3(refl3)
    flagged3 = refl3.get_flags(refl3.flags.centroid_outlier)
    assert (flagged3 & truth3).count(True) == truth3.count(True)
