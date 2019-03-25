# -*- coding: utf-8 -*-
"""
This module defines a number of general plots, which may be relevant to
for reports of several programs.
"""
from scitbx.array_family import flex
from cctbx import uctbx

def scale_rmerge_vs_batch_plot(batch_manager, rmerge_vs_b, scales_vs_b=None):
    reduced_batches = batch_manager.reduced_batches
    shapes, annotations, text = batch_manager.batch_plot_shapes_and_annotations()
    if len(annotations) > 30:
        # at a certain point the annotations become unreadable
        annotations = None

    return {
        "scale_rmerge_vs_batch": {
            "data": [
                (
                    {
                        "x": reduced_batches,
                        "y": scales_vs_b,
                        "type": "scatter",
                        "name": "Scale",
                        "opacity": 0.75,
                        "text": text,
                    }
                    if scales_vs_b is not None
                    else {}
                ),
                {
                    "x": reduced_batches,
                    "y": rmerge_vs_b,
                    "yaxis": "y2",
                    "type": "scatter",
                    "name": "Rmerge",
                    "opacity": 0.75,
                    "text": text,
                },
            ],
            "layout": {
                "title": "Scale and Rmerge vs batch",
                "xaxis": {"title": "N"},
                "yaxis": {"title": "Scale", "rangemode": "tozero"},
                "yaxis2": {
                    "title": "Rmerge",
                    "overlaying": "y",
                    "side": "right",
                    "rangemode": "tozero",
                },
                "shapes": shapes,
                "annotations": annotations,
            },
        }
    }


def i_over_sig_i_vs_batch_plot(batch_manager, i_sig_i_vs_batch):

    reduced_batches = batch_manager.reduced_batches
    shapes, annotations, _ = batch_manager.batch_plot_shapes_and_annotations()
    if len(annotations) > 30:
        # at a certain point the annotations become unreadable
        annotations = None

    return {
        "i_over_sig_i_vs_batch": {
            "data": [
                {
                    "x": reduced_batches,
                    "y": i_sig_i_vs_batch,
                    "type": "scatter",
                    "name": "I/sigI vs batch",
                    "opacity": 0.75,
                }
            ],
            "layout": {
                "title": "<I/sig(I)> vs batch",
                "xaxis": {"title": "N"},
                "yaxis": {"title": "<I/sig(I)>", "rangemode": "tozero"},
                "shapes": shapes,
                "annotations": annotations,
            },
        }
    }

def cc_one_half_plot(dataset_statistics, method=None, is_centric=False):

    if method == "sigma_tau":
        cc_one_half_bins = [
            bin_stats.cc_one_half_sigma_tau for bin_stats in dataset_statistics.bins
        ]
        cc_one_half_critical_value_bins = [
            bin_stats.cc_one_half_sigma_tau_critical_value
            for bin_stats in dataset_statistics.bins
        ]
    else:
        cc_one_half_bins = [
            bin_stats.cc_one_half for bin_stats in dataset_statistics.bins
        ]
        cc_one_half_critical_value_bins = [
            bin_stats.cc_one_half_critical_value
            for bin_stats in dataset_statistics.bins
        ]
    cc_anom_bins = [bin_stats.cc_anom for bin_stats in dataset_statistics.bins]
    cc_anom_critical_value_bins = [
        bin_stats.cc_anom_critical_value for bin_stats in dataset_statistics.bins
    ]

    d_star_sq_bins = [
        (1 / bin_stats.d_min ** 2) for bin_stats in dataset_statistics.bins
    ]
    d_star_sq_tickvals, d_star_sq_ticktext = _d_star_sq_to_d_ticks(
        d_star_sq_bins, nticks=5
    )

    return {
        "cc_one_half": {
            "data": [
                {
                    "x": d_star_sq_bins,  # d_star_sq
                    "y": cc_one_half_bins,
                    "type": "scatter",
                    "name": "CC-half",
                    "mode": "lines",
                    "line": {"color": "rgb(31, 119, 180)"},
                },
                {
                    "x": d_star_sq_bins,  # d_star_sq
                    "y": cc_one_half_critical_value_bins,
                    "type": "scatter",
                    "name": "CC-half critical value (p=0.01)",
                    "line": {"color": "rgb(31, 119, 180)", "dash": "dot"},
                },
                (
                    {
                        "x": d_star_sq_bins,  # d_star_sq
                        "y": cc_anom_bins,
                        "type": "scatter",
                        "name": "CC-anom",
                        "mode": "lines",
                        "line": {"color": "rgb(255, 127, 14)"},
                    }
                    if not is_centric
                    else {}
                ),
                (
                    {
                        "x": d_star_sq_bins,  # d_star_sq
                        "y": cc_anom_critical_value_bins,
                        "type": "scatter",
                        "name": "CC-anom critical value (p=0.01)",
                        "mode": "lines",
                        "line": {"color": "rgb(255, 127, 14)", "dash": "dot"},
                    }
                    if not is_centric
                    else {}
                ),
            ],
            "layout": {
                "title": "CC-half vs resolution",
                "xaxis": {
                    "title": u"Resolution (Å)",
                    "tickvals": d_star_sq_tickvals,
                    "ticktext": d_star_sq_ticktext,
                },
                "yaxis": {
                    "title": "CC-half",
                    "range": [min(cc_one_half_bins + cc_anom_bins + [0]), 1],
                },
            },
            "help": """\
The correlation coefficients, CC1/2, between random half-datasets. A correlation
coefficient of +1 indicates good correlation, and 0 indicates no correlation.
CC1/2 is typically close to 1 at low resolution, falling off to close to zero at
higher resolution. A typical resolution cutoff based on CC1/2 is around 0.3-0.5.

[1] Karplus, P. A., & Diederichs, K. (2012). Science, 336(6084), 1030-1033.
    https://doi.org/10.1126/science.1218231
[2] Diederichs, K., & Karplus, P. A. (2013). Acta Cryst D, 69(7), 1215-1222.
    https://doi.org/10.1107/S0907444913001121
[3] Evans, P. R., & Murshudov, G. N. (2013). Acta Cryst D, 69(7), 1204-1214.
    https://doi.org/10.1107/S0907444913000061
""",
        }
    }


def statistics_tables(dataset_statistics):
    result = dataset_statistics
    resolution_binned_table = [
        (
            "d_max",
            "d_min",
            "n_obs",
            "n_uniq",
            "mult",
            "comp",
            "&ltI&gt",
            "&ltI/sI&gt",
            "r_merge",
            "r_meas",
            "r_pim",
            "cc1/2",
            "cc_anom",
        )
    ]
    for bin_stats in result.bins:
        resolution_binned_table.append(tuple(bin_stats.format().split()))
    result = result.overall
    summary_table = [
        ("Resolution", "{0:.3f} - {1:.3f}".format(result.d_max, result.d_min)),
        ("Observations", result.n_obs),
        ("Unique Reflections", result.n_uniq),
        ("Redundancy", "{:.2f}".format(result.mean_redundancy)),
        ("Completeness", "{:.2f}".format(result.completeness * 100)),
        ("Mean intensity", "{:.1f}".format(result.i_mean)),
        ("Mean I/sigma(I)", "{:.1f}".format(result.i_over_sigma_mean)),
        ("R-merge", "{:.4f}".format(result.r_merge)),
        ("R-meas", "{:.4f}".format(result.r_meas)),
        ("R-pim", "{:.4f}".format(result.r_pim)),
    ]
    return (summary_table, resolution_binned_table)

def _d_star_sq_to_d_ticks(d_star_sq, nticks):
    d_spacings = uctbx.d_star_sq_as_d(flex.double(d_star_sq))
    min_d_star_sq = min(d_star_sq)
    dstep = (max(d_star_sq) - min_d_star_sq) / nticks
    tickvals = list(min_d_star_sq + (i * dstep) for i in range(nticks))
    ticktext = ["%.2f" % (uctbx.d_star_sq_as_d(dsq)) for dsq in tickvals]
    return tickvals, ticktext
