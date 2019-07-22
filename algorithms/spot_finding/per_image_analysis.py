from __future__ import absolute_import, division, print_function

import math
import sys

from cctbx import sgtbx, uctbx
from dials.algorithms.integration import filtering
from dials.array_family import flex
from libtbx import group_args, table_utils
from libtbx.math_utils import nearest_integer as nint
from scitbx import matrix


class slot(object):
    def __init__(self, d_min, d_max):
        self.d_min = d_min
        self.d_max = d_max


class binner_equal_population(object):
    def __init__(self, d_star_sq, target_n_per_bin=20, max_slots=20, min_slots=5):
        n_slots = len(d_star_sq) // target_n_per_bin
        if max_slots is not None:
            n_slots = min(n_slots, max_slots)
        if min_slots is not None:
            n_slots = max(n_slots, min_slots)
        self.bins = []
        n_per_bin = len(d_star_sq) / n_slots
        d_star_sq_sorted = flex.sorted(d_star_sq)
        d_sorted = uctbx.d_star_sq_as_d(d_star_sq_sorted)
        d_max = d_sorted[0]
        for i in range(n_slots):
            d_min = d_sorted[nint((i + 1) * n_per_bin) - 1]
            self.bins.append(slot(d_min, d_max))
            d_max = d_min
            # print self.bins[-1].d_max, self.bins[-1].d_min


class binner_d_star_cubed(object):
    def __init__(self, d_spacings, target_n_per_bin=25, max_slots=40, min_slots=20):
        d_spacings = flex.double(list(set(d_spacings)))
        d_spacings_sorted = flex.sorted(d_spacings, reverse=True)
        d_star_cubed_sorted = flex.pow(1 / d_spacings_sorted, 3)

        # choose bin volume such that lowest resolution shell contains 5% of the
        # spots, or 25, whichever is greater
        low_res_count = int(
            math.ceil(
                min(
                    max(target_n_per_bin, 0.05 * len(d_spacings)),
                    0.25 * len(d_spacings),
                )
            )
        )
        bin_step = d_star_cubed_sorted[low_res_count] - d_star_cubed_sorted[0]
        assert bin_step > 0
        n_slots = int(
            math.ceil((d_star_cubed_sorted[-1] - d_star_cubed_sorted[0]) / bin_step)
        )

        # n_slots = len(d_spacings_sorted)//target_n_per_bin
        if max_slots is not None:
            n_slots = min(n_slots, max_slots)
        if min_slots is not None:
            n_slots = max(n_slots, min_slots)
        bin_step = (d_star_cubed_sorted[-1] - d_star_cubed_sorted[0]) / n_slots

        self.bins = []
        ds3_max = d_star_cubed_sorted[0]
        for i in range(n_slots):
            ds3_min = d_star_cubed_sorted[0] + (i + 1) * bin_step
            self.bins.append(slot(1 / ds3_min ** (1 / 3), 1 / ds3_max ** (1 / 3)))
            ds3_max = ds3_min


def map_to_reciprocal_space(reflections, imageset):
    detector = imageset.get_detector()
    scan = imageset.get_scan()
    goniometer = imageset.get_goniometer()
    beam = imageset.get_beam()

    reflections.centroid_px_to_mm(detector, scan)
    reflections.map_centroids_to_reciprocal_space(detector, beam, goniometer)

    return reflections


def get_histogram(d_star_sq, target_n_per_bin=20, max_slots=20, min_slots=5):
    n_slots = len(d_star_sq) // target_n_per_bin
    n_slots = min(n_slots, max_slots)
    n_slots = max(n_slots, min_slots)
    return flex.histogram(d_star_sq, n_slots=n_slots)


def outlier_rejection(reflections):
    # http://scripts.iucr.org/cgi-bin/paper?ba0032
    if len(reflections) == 1:
        return reflections
    intensities = reflections["intensity.sum.value"]
    variances = reflections["intensity.sum.variance"]

    i_max = flex.max_index(intensities)

    sel = flex.bool(len(reflections), True)
    sel[i_max] = False

    i_test = intensities[i_max]
    var_test = variances[i_max]

    intensities_subset = intensities.select(sel)
    var_subset = variances.select(sel)

    var_prior = var_test + 1 / flex.sum(1 / var_subset)
    p_prior = (
        1
        / math.sqrt(2 * math.pi * var_prior)
        * math.exp(-(i_test - flex.mean(intensities_subset)) ** 2 / (2 * var_prior))
    )
    # print p_prior

    if p_prior > 1e-10:
        return reflections

    return outlier_rejection(reflections.select(sel))


def wilson_outliers(reflections, ice_sel=None, p_cutoff=1e-2):
    # http://scripts.iucr.org/cgi-bin/paper?ba0032
    if ice_sel is None:
        ice_sel = flex.bool(len(reflections), False)

    E_cutoff = math.sqrt(-math.log(p_cutoff))
    intensities = reflections["intensity.sum.value"]

    Sigma_n = flex.mean(intensities.select(~ice_sel))
    normalised_amplitudes = flex.sqrt(intensities) / math.sqrt(Sigma_n)

    outliers = normalised_amplitudes >= E_cutoff

    if outliers.count(True):
        # iterative outlier rejection
        inliers = ~outliers
        outliers.set_selected(
            inliers,
            wilson_outliers(reflections.select(inliers), ice_sel.select(inliers)),
        )

    return outliers


def estimate_resolution_limit(reflections, imageset, ice_sel=None, plot_filename=None):
    if ice_sel is None:
        ice_sel = flex.bool(len(reflections), False)

    d_star_sq = flex.pow2(reflections["rlp"].norms())
    d_spacings = uctbx.d_star_sq_as_d(d_star_sq)

    intensities = reflections["intensity.sum.value"]
    variances = reflections["intensity.sum.variance"]

    sel = variances > 0
    intensities = intensities.select(sel)
    variances = variances.select(sel)
    ice_sel = ice_sel.select(sel)

    i_over_sigi = intensities / flex.sqrt(variances)
    log_i_over_sigi = flex.log(i_over_sigi)

    fit = flex.linear_regression(
        d_star_sq.select(~ice_sel), log_i_over_sigi.select(~ice_sel)
    )
    m = fit.slope()
    c = fit.y_intercept()

    log_i_sigi_lower = flex.double()
    d_star_sq_lower = flex.double()
    log_i_sigi_upper = flex.double()
    d_star_sq_upper = flex.double()

    binner = binner_equal_population(
        d_star_sq, target_n_per_bin=20, max_slots=20, min_slots=5
    )

    outliers_all = flex.bool(len(reflections), False)

    low_percentile_limit = 0.1
    upper_percentile_limit = 1 - low_percentile_limit
    for i_slot, slot in enumerate(binner.bins):
        sel_all = (d_spacings < slot.d_max) & (d_spacings >= slot.d_min)
        sel = ~(ice_sel) & sel_all

        if sel.count(True) == 0:
            continue

        outliers = wilson_outliers(
            reflections.select(sel_all), ice_sel=ice_sel.select(sel_all)
        )
        outliers_all.set_selected(sel_all, outliers)

        isel = sel_all.iselection().select(~(outliers) & ~(ice_sel).select(sel_all))
        log_i_over_sigi_sel = log_i_over_sigi.select(isel)
        d_star_sq_sel = d_star_sq.select(isel)

        perm = flex.sort_permutation(log_i_over_sigi_sel)
        i_lower = perm[int(math.floor(low_percentile_limit * len(perm)))]
        i_upper = perm[int(math.floor(upper_percentile_limit * len(perm)))]
        log_i_sigi_lower.append(log_i_over_sigi_sel[i_lower])
        log_i_sigi_upper.append(log_i_over_sigi_sel[i_upper])
        d_star_sq_upper.append(d_star_sq_sel[i_lower])
        d_star_sq_lower.append(d_star_sq_sel[i_upper])

    fit_upper = flex.linear_regression(d_star_sq_upper, log_i_sigi_upper)
    m_upper = fit_upper.slope()
    c_upper = fit_upper.y_intercept()
    fit_lower = flex.linear_regression(d_star_sq_lower, log_i_sigi_lower)
    m_lower = fit_lower.slope()
    c_lower = fit_lower.y_intercept()

    # fit_upper.show_summary()
    # fit_lower.show_summary()

    if m_upper == m_lower:
        intersection = (-1, -1)
        resolution_estimate = -1
        inside = flex.bool(len(d_star_sq), False)

    else:
        # http://en.wikipedia.org/wiki/Line%E2%80%93line_intersection#Given_the_equations_of_the_lines
        # with:
        # a_ = m_upper
        # b_ = m_lower
        # c_ = c_upper
        # d_ = c_lower
        # intersection == ((d_ - c_) / (a_ - b_), (a_ * d_ - b_ * c_) / (a_ - b_))
        intersection = (
            (c_lower - c_upper) / (m_upper - m_lower),
            (m_upper * c_lower - m_lower * c_upper) / (m_upper - m_lower),
        )

        # inside = points_inside_envelope(
        # d_star_sq, log_i_over_sigi, m_upper, c_upper, m_lower, c_lower)

        inside = points_below_line(d_star_sq, log_i_over_sigi, m_upper, c_upper)
        inside = inside & ~outliers_all

        if inside.count(True) > 0:
            d_star_sq_estimate = flex.max(d_star_sq.select(inside))
            # d_star_sq_estimate = intersection[0]
            resolution_estimate = uctbx.d_star_sq_as_d(d_star_sq_estimate)
        else:
            resolution_estimate = -1

    # resolution_estimate = max(resolution_estimate, flex.min(d_spacings))

    if plot_filename is not None:
        from matplotlib import pyplot

        fig = pyplot.figure()
        ax = fig.add_subplot(1, 1, 1)
        ax.scatter(d_star_sq, log_i_over_sigi, marker="+")
        ax.scatter(
            d_star_sq.select(inside),
            log_i_over_sigi.select(inside),
            marker="+",
            color="green",
        )
        ax.scatter(
            d_star_sq.select(ice_sel),
            log_i_over_sigi.select(ice_sel),
            marker="+",
            color="black",
        )
        ax.scatter(
            d_star_sq.select(outliers_all),
            log_i_over_sigi.select(outliers_all),
            marker="+",
            color="grey",
        )
        ax.scatter(d_star_sq_upper, log_i_sigi_upper, marker="+", color="red")
        ax.scatter(d_star_sq_lower, log_i_sigi_lower, marker="+", color="red")

        if intersection[0] <= ax.get_xlim()[1] and intersection[1] <= ax.get_ylim()[1]:
            ax.scatter(
                [intersection[0]], [intersection[1]], marker="x", s=50, color="b"
            )
        # ax.hexbin(d_star_sq, log_i_over_sigi, gridsize=30)
        xlim = pyplot.xlim()
        ax.plot(xlim, [(m * x + c) for x in xlim])
        ax.plot(xlim, [(m_upper * x + c_upper) for x in xlim], color="red")
        ax.plot(xlim, [(m_lower * x + c_lower) for x in xlim], color="red")
        ax.set_xlabel("d_star_sq")
        ax.set_ylabel("ln(I/sigI)")
        ax.set_xlim((max(-xlim[1], -0.05), xlim[1]))
        ax.set_ylim((0, ax.get_ylim()[1]))

        for i_slot, slot in enumerate(binner.bins):
            if i_slot == 0:
                ax.vlines(
                    uctbx.d_as_d_star_sq(slot.d_max),
                    0,
                    ax.get_ylim()[1],
                    linestyle="dotted",
                    color="grey",
                )
            ax.vlines(
                uctbx.d_as_d_star_sq(slot.d_min),
                0,
                ax.get_ylim()[1],
                linestyle="dotted",
                color="grey",
            )

        ax_ = ax.twiny()  # ax2 is responsible for "top" axis and "right" axis
        xticks = ax.get_xticks()
        xticks_d = [uctbx.d_star_sq_as_d(ds2) if ds2 > 0 else 0 for ds2 in xticks]
        ax_.set_xticks(xticks)
        ax_.set_xlim(ax.get_xlim())
        ax_.set_xlabel(r"Resolution ($\AA$)")
        ax_.set_xticklabels(["%.1f" % d for d in xticks_d])
        pyplot.savefig(plot_filename)
        pyplot.close()

    return resolution_estimate


def estimate_resolution_limit_distl_method1(reflections, imageset, plot_filename=None):
    # Implementation of Method 1 (section 2.4.4) of:
    # Z. Zhang, N. K. Sauter, H. van den Bedem, G. Snell and A. M. Deacon
    # J. Appl. Cryst. (2006). 39, 112-119
    # https://doi.org/10.1107/S0021889805040677

    variances = reflections["intensity.sum.variance"]

    sel = variances > 0
    reflections = reflections.select(sel)
    d_star_sq = flex.pow2(reflections["rlp"].norms())
    d_spacings = uctbx.d_star_sq_as_d(d_star_sq)
    d_star_cubed = flex.pow(reflections["rlp"].norms(), 3)

    step = 2
    while len(reflections) / step > 40:
        step += 1

    order = flex.sort_permutation(d_spacings, reverse=True)

    ds3_subset = flex.double()
    d_subset = flex.double()
    for i in range(len(reflections) // step):
        ds3_subset.append(d_star_cubed[order[i * step]])
        d_subset.append(d_spacings[order[i * step]])

    x = flex.double(range(len(ds3_subset)))

    # (i)
    # Usually, Pm is the last point, that is, m = n. But m could be smaller than
    # n if an unusually high number of spots are detected around a certain
    # intermediate resolution. In that case, our search for the image resolution
    # does not go outside the spot 'bump;. This is particularly useful when
    # ice-rings are present.

    slopes = (ds3_subset[1:] - ds3_subset[0]) / (x[1:] - x[0])
    skip_first = 3
    p_m = flex.max_index(slopes[skip_first:]) + 1 + skip_first

    # (ii)

    x1 = matrix.col((0, ds3_subset[0]))
    x2 = matrix.col((p_m, ds3_subset[p_m]))

    gaps = flex.double([0])
    v = matrix.col(((x2[1] - x1[1]), -(x2[0] - x1[0]))).normalize()

    for i in range(1, p_m):
        x0 = matrix.col((i, ds3_subset[i]))
        r = x1 - x0
        g = abs(v.dot(r))
        gaps.append(g)

    mv = flex.mean_and_variance(gaps)
    s = mv.unweighted_sample_standard_deviation()

    # (iii)

    p_k = flex.max_index(gaps)
    g_k = gaps[p_k]
    p_g = p_k
    for i in range(p_k + 1, len(gaps)):
        g_i = gaps[i]
        if g_i > (g_k - 0.5 * s):
            p_g = i

    d_g = d_subset[p_g]

    noisiness = 0
    n = len(ds3_subset)
    for i in range(n - 1):
        for j in range(i + 1, n - 1):
            if slopes[i] >= slopes[j]:
                noisiness += 1
    noisiness /= (n - 1) * (n - 2) / 2

    if plot_filename is not None:
        from matplotlib import pyplot

        fig = pyplot.figure()
        ax = fig.add_subplot(1, 1, 1)
        ax.scatter(range(len(ds3_subset)), ds3_subset)
        ax.set_ylabel("D^-3")
        xlim = pyplot.xlim()
        ylim = pyplot.ylim()
        ax.vlines(p_g, ylim[0], ylim[1], colors="red")
        pyplot.xlim(0, xlim[1])
        pyplot.ylim(0, ylim[1])
        pyplot.savefig(plot_filename)
        pyplot.close()

    return d_g, noisiness


def estimate_resolution_limit_distl_method2(reflections, imageset, plot_filename=None):
    # Implementation of Method 2 (section 2.4.4) of:
    # Z. Zhang, N. K. Sauter, H. van den Bedem, G. Snell and A. M. Deacon
    # J. Appl. Cryst. (2006). 39, 112-119
    # https://doi.org/10.1107/S0021889805040677

    variances = reflections["intensity.sum.variance"]

    sel = variances > 0
    reflections = reflections.select(sel)
    d_star_sq = flex.pow2(reflections["rlp"].norms())
    d_spacings = uctbx.d_star_sq_as_d(d_star_sq)

    binner = binner_d_star_cubed(d_spacings)

    bin_counts = flex.size_t()

    for i_slot, slot in enumerate(binner.bins):
        sel_all = (d_spacings < slot.d_max) & (d_spacings >= slot.d_min)
        sel = sel_all

        bin_counts.append(sel.count(True))

    # print list(bin_counts)
    t0 = (bin_counts[0] + bin_counts[1]) / 2

    mu = 0.15

    for i in range(len(bin_counts) - 1):
        tj = bin_counts[i]
        tj1 = bin_counts[i + 1]
        if (tj < (mu * t0)) and (tj1 < (mu * t0)):
            break

    d_min = binner.bins[i].d_min
    noisiness = 0
    m = len(bin_counts)
    for i in range(m):
        for j in range(i + 1, m):
            if bin_counts[i] <= bin_counts[j]:
                noisiness += 1
    noisiness /= 0.5 * m * (m - 1)

    if plot_filename is not None:
        from matplotlib import pyplot

        fig = pyplot.figure()
        ax = fig.add_subplot(1, 1, 1)
        ax.scatter(range(len(bin_counts)), bin_counts)
        ax.set_ylabel("number of spots in shell")
        xlim = pyplot.xlim()
        ylim = pyplot.ylim()
        ax.vlines(i, ylim[0], ylim[1], colors="red")
        pyplot.xlim(0, xlim[1])
        pyplot.ylim(0, ylim[1])
        pyplot.savefig(plot_filename)
        pyplot.close()

    return d_min, noisiness


def points_below_line(d_star_sq, log_i_over_sigi, m, c):
    p1 = matrix.col((0, c))
    p2 = matrix.col((1, m * 1 + c))

    def side(p1, p2, p):
        diff = p2 - p1
        perp = matrix.col((-diff[1], diff[0]))
        # print p, p1, p2, perp
        d = (p - p1).dot(perp)
        # print d
        return math.copysign(1, d)

    inside = flex.bool(len(d_star_sq), False)
    for i, (x, y) in enumerate(zip(d_star_sq, log_i_over_sigi)):
        p = matrix.col((x, y))
        if side(p1, p2, p) < 0:
            inside[i] = True

    return inside


def points_inside_envelope(
    d_star_sq, log_i_over_sigi, m_upper, c_upper, m_lower, c_lower
):

    return points_below_line(
        d_star_sq, log_i_over_sigi, m_upper, c_upper
    ) & ~points_below_line(d_star_sq, log_i_over_sigi, m_lower, c_lower)


def ice_rings_selection(reflections, width=0.004):
    d_star_sq = flex.pow2(reflections["rlp"].norms())
    d_spacings = uctbx.d_star_sq_as_d(d_star_sq)

    unit_cell = uctbx.unit_cell((4.498, 4.498, 7.338, 90, 90, 120))
    space_group = sgtbx.space_group_info(number=194).group()

    if d_spacings:
        ice_filter = filtering.PowderRingFilter(
            unit_cell, space_group, flex.min(d_spacings), width
        )

        ice_sel = ice_filter(d_spacings)

        return ice_sel
    else:
        return None


def resolution_histogram(reflections, imageset, plot_filename=None):
    d_star_sq = flex.pow2(reflections["rlp"].norms())
    hist = get_histogram(d_star_sq)

    if plot_filename is not None:
        from matplotlib import pyplot

        fig = pyplot.figure()
        ax = fig.add_subplot(1, 1, 1)
        ax.bar(
            hist.slot_centers() - 0.5 * hist.slot_width(),
            hist.slots(),
            width=hist.slot_width(),
        )
        ax.set_xlabel("d_star_sq")
        ax.set_ylabel("Frequency")

        ax_ = ax.twiny()  # ax2 is responsible for "top" axis and "right" axis
        xticks = ax.get_xticks()
        xticks_d = [uctbx.d_star_sq_as_d(ds2) if ds2 > 0 else 0 for ds2 in xticks]
        ax_.set_xticks(xticks)
        ax_.set_xlim(ax.get_xlim())
        ax_.set_xlabel(r"Resolution ($\AA$)")
        ax_.set_xticklabels(["%.1f" % d for d in xticks_d])
        pyplot.savefig(plot_filename)
        pyplot.close()


def log_sum_i_sigi_vs_resolution(reflections, imageset, plot_filename=None):
    d_star_sq = flex.pow2(reflections["rlp"].norms())
    hist = get_histogram(d_star_sq)

    intensities = reflections["intensity.sum.value"]
    variances = reflections["intensity.sum.variance"]

    sel = variances > 0
    intensities = intensities.select(sel)
    variances = intensities.select(sel)

    i_over_sigi = intensities / flex.sqrt(variances)
    # log_i_over_sigi = flex.log(i_over_sigi)

    slots = []
    for slot in hist.slot_infos():
        sel = (d_star_sq > slot.low_cutoff) & (d_star_sq < slot.high_cutoff)
        if sel.count(True) > 0:
            slots.append(math.log(flex.sum(i_over_sigi.select(sel))))
        else:
            slots.append(0)

    if plot_filename is not None:
        from matplotlib import pyplot

        fig = pyplot.figure()
        ax = fig.add_subplot(1, 1, 1)
        # ax.bar(hist.slot_centers()-0.5*hist.slot_width(), hist.slots(),
        ax.scatter(
            hist.slot_centers() - 0.5 * hist.slot_width(),
            slots,
            s=20,
            color="blue",
            marker="o",
            alpha=0.5,
        )
        ax.set_xlabel("d_star_sq")
        ax.set_ylabel("ln(sum(I/sigI))")

        ax_ = ax.twiny()  # ax2 is responsible for "top" axis and "right" axis
        xticks = ax.get_xticks()
        xticks_d = [uctbx.d_star_sq_as_d(ds2) if ds2 > 0 else 0 for ds2 in xticks]
        ax_.set_xticks(xticks)
        ax_.set_xlim(ax.get_xlim())
        ax_.set_xlabel(r"Resolution ($\AA$)")
        ax_.set_xticklabels(["%.1f" % d for d in xticks_d])
        pyplot.savefig(plot_filename)
        pyplot.close()


def plot_ordered_d_star_sq(reflections, imageset):
    from matplotlib import pyplot

    d_star_sq = flex.pow2(reflections["rlp"].norms())

    perm = flex.sort_permutation(d_star_sq)
    pyplot.scatter(list(range(len(perm))), list(d_star_sq.select(perm)), marker="+")
    pyplot.show()


def stats_single_image(
    imageset,
    reflections,
    i=None,
    resolution_analysis=True,
    plot=False,
    filter_ice=True,
    ice_rings_width=0.004,
):
    reflections = map_to_reciprocal_space(reflections, imageset)
    if plot and i is not None:
        filename = "i_over_sigi_vs_resolution_%d.png" % (i + 1)
        hist_filename = "spot_count_vs_resolution_%d.png" % (i + 1)
        extra_filename = "log_sum_i_sigi_vs_resolution_%d.png" % (i + 1)
        distl_method_1_filename = "distl_method_1_%d.png" % (i + 1)
        distl_method_2_filename = "distl_method_2_%d.png" % (i + 1)
    else:
        filename = None
        hist_filename = None
        extra_filename = None
        distl_method_1_filename = None
        distl_method_2_filename = None

    d_star_sq = flex.pow2(reflections["rlp"].norms())
    d_spacings = uctbx.d_star_sq_as_d(d_star_sq)

    # plot_ordered_d_star_sq(reflections, imageset)
    reflections_all = reflections
    reflections_no_ice = reflections_all
    ice_sel = None
    if filter_ice:
        ice_sel = ice_rings_selection(reflections_all, width=ice_rings_width)
        if ice_sel is not None:
            reflections_no_ice = reflections_all.select(~ice_sel)
    n_spots_total = len(reflections_all)
    n_spots_no_ice = len(reflections_no_ice)
    n_spot_4A = (d_spacings > 4).count(True)
    intensities = reflections_no_ice["intensity.sum.value"]
    total_intensity = flex.sum(intensities)
    # print i
    if hist_filename is not None:
        resolution_histogram(reflections, imageset, plot_filename=hist_filename)
    if extra_filename is not None:
        log_sum_i_sigi_vs_resolution(
            reflections, imageset, plot_filename=extra_filename
        )
    if resolution_analysis and n_spots_no_ice > 10:
        estimated_d_min = estimate_resolution_limit(
            reflections_all, imageset, ice_sel=ice_sel, plot_filename=filename
        )
        d_min_distl_method_1, noisiness_method_1 = estimate_resolution_limit_distl_method1(
            reflections_all, imageset, plot_filename=distl_method_1_filename
        )
        d_min_distl_method_2, noisiness_method_2 = estimate_resolution_limit_distl_method2(
            reflections_all, imageset, plot_filename=distl_method_2_filename
        )
    else:
        estimated_d_min = -1.0
        d_min_distl_method_1 = -1.0
        noisiness_method_1 = -1.0
        d_min_distl_method_2 = -1.0
        noisiness_method_2 = -1.0

    return group_args(
        n_spots_total=n_spots_total,
        n_spots_no_ice=n_spots_no_ice,
        n_spots_4A=n_spot_4A,
        total_intensity=total_intensity,
        estimated_d_min=estimated_d_min,
        d_min_distl_method_1=d_min_distl_method_1,
        noisiness_method_1=noisiness_method_1,
        d_min_distl_method_2=d_min_distl_method_2,
        noisiness_method_2=noisiness_method_2,
    )


def stats_imageset(imageset, reflections, resolution_analysis=True, plot=False):
    n_spots_total = []
    n_spots_no_ice = []
    n_spots_4A = []
    total_intensity = []
    estimated_d_min = []
    d_min_distl_method_1 = []
    d_min_distl_method_2 = []
    noisiness_method_1 = []
    noisiness_method_2 = []

    image_number = reflections["xyzobs.px.value"].parts()[2]
    image_number = flex.floor(image_number)

    try:
        start, end = imageset.get_array_range()
    except AttributeError:
        start = 0
    for i in range(len(imageset)):
        stats = stats_single_image(
            imageset[i : i + 1],
            reflections.select(image_number == i + start),
            i=i + start,
            resolution_analysis=resolution_analysis,
            plot=plot,
        )
        n_spots_total.append(stats.n_spots_total)
        n_spots_no_ice.append(stats.n_spots_no_ice)
        n_spots_4A.append(stats.n_spots_4A)
        total_intensity.append(stats.total_intensity)
        estimated_d_min.append(stats.estimated_d_min)
        d_min_distl_method_1.append(stats.d_min_distl_method_1)
        noisiness_method_1.append(stats.noisiness_method_1)
        d_min_distl_method_2.append(stats.d_min_distl_method_2)
        noisiness_method_2.append(stats.noisiness_method_2)

    return group_args(
        n_spots_total=n_spots_total,
        n_spots_no_ice=n_spots_no_ice,
        n_spots_4A=n_spots_4A,
        total_intensity=total_intensity,
        estimated_d_min=estimated_d_min,
        d_min_distl_method_1=d_min_distl_method_1,
        noisiness_method_1=noisiness_method_1,
        d_min_distl_method_2=d_min_distl_method_2,
        noisiness_method_2=noisiness_method_2,
    )


def table(stats, perm=None, n_rows=None):
    n_spots_total = stats.n_spots_total
    n_spots_no_ice = stats.n_spots_no_ice
    total_intensity = stats.total_intensity
    estimated_d_min = stats.estimated_d_min
    d_min_distl_method_1 = stats.d_min_distl_method_1
    d_min_distl_method_2 = stats.d_min_distl_method_2
    noisiness_method_1 = stats.noisiness_method_1
    noisiness_method_2 = stats.noisiness_method_2
    if hasattr(stats, "image"):
        image = stats.image
    else:
        image = flex.int(range(1, len(n_spots_total) + 1)).as_string()
    n_indexed = None
    fraction_indexed = None
    if hasattr(stats, "n_indexed"):
        n_indexed = stats.n_indexed
    if hasattr(stats, "fraction_indexed"):
        fraction_indexed = stats.fraction_indexed
    if flex.double(estimated_d_min).all_eq(-1):
        estimated_d_min = None
    if flex.double(d_min_distl_method_1).all_eq(-1):
        d_min_distl_method_1 = None
    if flex.double(d_min_distl_method_2).all_eq(-1):
        d_min_distl_method_2 = None

    rows = [["image", "#spots", "#spots_no_ice", "total_intensity"]]
    if estimated_d_min is not None:
        rows[0].append("d_min")
    if d_min_distl_method_1 is not None:
        rows[0].append("d_min (distl method 1)")
    if d_min_distl_method_2 is not None:
        rows[0].append("d_min (distl method 2)")
    if n_indexed is not None:
        rows[0].append("#indexed")
    if fraction_indexed is not None:
        rows[0].append("fraction_indexed")
    if perm is None:
        perm = list(range(len(n_spots_total)))
    if n_rows is not None:
        n_rows = min(n_rows, len(perm))
        perm = perm[:n_rows]
    for i_image in perm:
        d_min_str = ""
        method1_str = ""
        method2_str = ""
        if estimated_d_min is not None and estimated_d_min[i_image] > 0:
            d_min_str = "%.2f" % estimated_d_min[i_image]
        if d_min_distl_method_1 is not None and d_min_distl_method_1[i_image] > 0:
            method1_str = "%.2f" % d_min_distl_method_1[i_image]
            if noisiness_method_1 is not None:
                method1_str += " (%.2f)" % noisiness_method_1[i_image]
        if d_min_distl_method_2 is not None and d_min_distl_method_2[i_image] > 0:
            method2_str = "%.2f" % d_min_distl_method_2[i_image]
            if noisiness_method_2 is not None:
                method2_str += " (%.2f)" % noisiness_method_2[i_image]
        row = [
            image[i_image],
            str(n_spots_total[i_image]),
            str(n_spots_no_ice[i_image]),
            "%.0f" % total_intensity[i_image],
        ]
        if estimated_d_min is not None:
            row.append(d_min_str)
        if d_min_distl_method_1 is not None:
            row.append(method1_str)
        if d_min_distl_method_2 is not None:
            row.append(method2_str)
        if n_indexed is not None:
            row.append("%i" % n_indexed[i_image])
        if fraction_indexed is not None:
            row.append("%.2f" % fraction_indexed[i_image])
        rows.append(row)
    return rows


def print_table(stats, perm=None, n_rows=None, out=None):
    if out is None:
        out = sys.stdout

    rows = table(stats, perm=perm, n_rows=n_rows)
    print(
        table_utils.format(rows, has_header=True, prefix="| ", postfix=" |"), file=out
    )


def plot_stats(stats, filename="per_image_analysis.png"):
    from matplotlib import pyplot

    n_spots_total = flex.int(stats.n_spots_total)
    n_spots_no_ice = stats.n_spots_no_ice
    n_spots_4A = stats.n_spots_4A
    estimated_d_min = flex.double(stats.estimated_d_min)
    d_min_distl_method_1 = flex.double(stats.d_min_distl_method_1)
    d_min_distl_method_2 = flex.double(stats.d_min_distl_method_2)
    total_intensity = stats.total_intensity

    i_image = flex.int(list(range(1, len(n_spots_total) + 1)))

    _, (ax1, ax2, ax3) = pyplot.subplots(nrows=3)
    ax1.scatter(
        list(i_image),
        list(n_spots_total),
        s=5,
        color="orange",
        marker="o",
        alpha=0.4,
        label="#spots (total)",
    )
    if n_spots_4A is not None:
        ax1.scatter(
            list(i_image),
            n_spots_4A,
            s=5,
            color="green",
            marker="o",
            alpha=0.4,
            label=u"#spots (to 4\u00c5)",
        )
    ax1.scatter(
        list(i_image),
        n_spots_no_ice,
        s=5,
        color="blue",
        marker="o",
        alpha=0.4,
        label="#spots (no ice)",
    )
    ax1.set_xlabel("Image #")
    ax1.set_ylabel("# spots")
    ax1.set_xlim((0.0, len(n_spots_total)))
    ax1.set_ylim(bottom=-0.2)
    ax1.legend(bbox_to_anchor=(1.05, 0.5), loc="center left", borderaxespad=0.0)

    sel = (estimated_d_min < 50.0) & (n_spots_total > 20) & (estimated_d_min > 0)
    ax2.scatter(
        list(i_image.select(sel)),
        list(estimated_d_min.select(sel)),
        s=10,
        color="red",
        marker="^",
        alpha=0.5,
        label="Estimated d_min",
    )
    ax2.scatter(
        list(i_image.select(sel)),
        list(d_min_distl_method_1.select(sel)),
        s=10,
        color="purple",
        marker="x",
        alpha=0.5,
        label="d_min (distl method 1)",
    )
    ax2.scatter(
        list(i_image.select(sel)),
        list(d_min_distl_method_2.select(sel)),
        s=10,
        color="grey",
        marker="+",
        alpha=0.5,
        label="d_min (distl method 2)",
    )
    ax2.set_ylabel(u"Resolution (\u00c5)")
    ax2.set_xlim((0, len(n_spots_total)))
    ylim = ax2.get_ylim()
    ylim = (math.floor(ylim[0]), math.ceil(ylim[1]))
    ax2.set_ylim(ylim)
    ax2.invert_yaxis()
    ax2.legend(bbox_to_anchor=(1.05, 0.5), loc="center left", borderaxespad=0.0)

    ax3.scatter(
        i_image,
        total_intensity,
        s=5,
        color="grey",
        marker="o",
        alpha=0.4,
        label="Total intensity",
    )
    ax3.set_ylabel("Total intensity")
    ax3.set_xlabel("Image #")
    ax3.legend(bbox_to_anchor=(1.05, 0.5), loc="center left", borderaxespad=0.0)

    pyplot.savefig(filename, dpi=600, bbox_inches="tight")
