# -*- coding: utf-8 -*-
"""Definitions of functions and classes for scaling and filtering algorithm."""
from __future__ import absolute_import, division, print_function

import math
from collections import OrderedDict

from dials.report.plots import ResolutionPlotterMixin
from libtbx import phil
from scitbx.array_family import flex

phil_scope = phil.parse(
    """
filtering {
    method = None deltacchalf
        .type = choice
        .help = "Choice of whether to do any filtering cycles, default None."
    deltacchalf {
        max_cycles = 6
            .type = int(value_min=1)
        max_percent_removed = 10
            .type = float
        min_completeness = None
            .type = float(value_min=0, value_max=100)
            .help = "Desired minimum completeness, as a percentage (0 - 100)."
        mode = *dataset image_group
            .type = choice
            .help = "Perform analysis on whole datasets or batch groups"
        group_size = 10
            .type = int(value_min=1)
            .help = "The number of images to group together when calculating delta"
                    "cchalf in image_group mode"
        stdcutoff = 4.0
            .type = float
            .help = "Datasets with a delta cc half below (mean - stdcutoff*std) are removed"
    }
    output {
        scale_and_filter_results = "scale_and_filter_results.json"
            .type = str
            .help = "Filename for output json of scale and filter results."
    }
}
"""
)

ordinal = lambda n: "%d%s" % (
    n,
    "tsnrhtdd"[(n / 10 % 10 != 1) * (n % 10 < 4) * n % 10 :: 4],
)


def log_cycle_results(results, scaling_script, filter_script):
    """Log results from the scripts for this cycle and add to the results."""
    cycle_results = {"merging_stats": scaling_script.merging_statistics_result}

    if not results.get_cycle_results():
        results.initial_n_reflections = scaling_script.scaled_miller_array.size()

    cycle_results["delta_cc_half_values"] = filter_script.results_summary[
        "per_dataset_delta_cc_half_values"
    ]["delta_cc_half_values"]
    cycle_results["mean_cc_half"] = filter_script.results_summary["mean_cc_half"]
    removal_summary = filter_script.results_summary["dataset_removal"]
    if removal_summary["mode"] == "image_group":
        cycle_results["image_ranges_removed"] = removal_summary["image_ranges_removed"]
    cycle_results["removed_datasets"] = removal_summary["experiments_fully_removed"]

    cycle_results["n_removed"] = filter_script.results_summary["dataset_removal"][
        "n_reflections_removed"
    ]

    n_removed = (
        sum(res["n_removed"] for res in results.get_cycle_results())
        + cycle_results["n_removed"]
    )
    percent_removed = n_removed / results.initial_n_reflections * 100
    cycle_results["cumul_percent_removed"] = percent_removed

    results.add_cycle_results(cycle_results)
    return results


class AnalysisResults(object):
    """Class to store results from scaling and filtering."""

    def __init__(self):
        self.termination_reason = None
        self.cycle_results = []
        self.initial_n_reflections = None
        self.initial_expids_and_image_ranges = None
        self.expids_and_image_ranges = None
        self.final_stats = None

    def add_cycle_results(self, results_dict):
        """Add the results dict from a scale and filter cycle."""
        merging_stats_dict = self._parse_merging_stats(results_dict["merging_stats"])
        results_dict["merging_stats"] = merging_stats_dict
        self.cycle_results.append(results_dict)

    @staticmethod
    def _parse_merging_stats(merging_stats_obj):
        merging_stats = {}
        overall = merging_stats_obj.overall
        merging_stats["ccs"] = [b.cc_one_half for b in merging_stats_obj.bins]
        merging_stats["rmerge"] = [b.r_merge for b in merging_stats_obj.bins]
        merging_stats["rpim"] = [b.r_pim for b in merging_stats_obj.bins]
        merging_stats["d_min"] = [b.d_min for b in merging_stats_obj.bins]
        merging_stats["overall"] = OrderedDict(
            [
                ("cc_one_half", overall.cc_one_half),
                ("r_merge", overall.r_merge),
                ("r_pim", overall.r_pim),
                ("i_over_sigma_mean", overall.i_over_sigma_mean),
                ("completeness", 100 * overall.completeness),
                ("n_obs", overall.n_obs),
            ]
        )
        return merging_stats

    def get_cycle_results(self):
        """Get the results from all cycles."""
        return self.cycle_results

    def get_last_cycle_results(self):
        """Get the results from the latest recorded cycle."""
        return self.cycle_results[-1]

    def add_final_stats(self, final_stats):
        """Add additional final merging stats from final rescale."""
        self.final_stats = self._parse_merging_stats(final_stats)

    def get_merging_stats(self):
        """Get all merging stats, including additional final stats if present."""
        stats = [res["merging_stats"] for res in self.cycle_results]
        if self.final_stats:
            stats += [self.final_stats]
        return stats

    def finish(self, termination_reason):
        """Set the termination reason/"""
        assert termination_reason in [
            "no_more_removed",
            "max_cycles",
            "max_percent_removed",
            "below_completeness_limit",
        ]
        self.termination_reason = termination_reason

    def to_dict(self):
        """Return the stored data as a dictionary."""
        return {
            "termination_reason": self.termination_reason,
            "initial_n_reflections": self.initial_n_reflections,
            "initial_expids_and_image_ranges": self.initial_expids_and_image_ranges,
            "cycle_results": OrderedDict(
                (i + 1, val) for i, val in enumerate(self.cycle_results)
            ),
            "final_stats": self.final_stats,
        }

    @staticmethod
    def from_dict(dictionary):
        """Configure the class from its dictionary form."""
        results = AnalysisResults()
        results.termination_reason = dictionary["termination_reason"]
        results.initial_expids_and_image_ranges = dictionary[
            "initial_expids_and_image_ranges"
        ]
        results.cycle_results = [
            dictionary["cycle_results"][key]
            for key in sorted(dictionary["cycle_results"].iterkeys())
        ]
        results.initial_n_reflections = dictionary["initial_n_reflections"]
        results.final_stats = dictionary["final_stats"]
        return results

    def make_summary(self):
        """Make summary of results."""
        msg = "\nSummary of data removed:\n"
        for i, res in enumerate(self.get_cycle_results()):
            msg += "Cycle number: %s\n" % (i + 1)
            if "image_ranges_removed" in res:
                if res["image_ranges_removed"]:
                    removed = "\n    ".join(
                        str(t[0]) + ", dataset " + str(t[1])
                        for t in res["image_ranges_removed"]
                    )
                    msg += "  Removed image ranges: \n    %s" % removed
            else:
                if res["removed_datasets"]:
                    msg += "  Removed datasets: %s\n" % res["removed_datasets"]
            msg += (
                "  cumulative %% of reflections removed: %.3f\n"
                % res["cumul_percent_removed"]
            )
        return msg


color_list = [
    "#F44336",
    "#FFC107",
    "#FFEB3B",
    "#8BC34A",
    "#03A9F4",
    "#3F51B5",
    "#607D8B",
]


def make_scaling_filtering_plots(data):
    """Make the plots for scaling and filtering."""
    d = make_filtering_merging_stats_plots(data["merging_stats"])
    d.update(make_histogram_plots(data["cycle_results"]))
    if data["mode"] == "image_group":
        d.update(
            make_reduction_plots(
                data["initial_expids_and_image_ranges"], data["expids_and_image_ranges"]
            )
        )
    return d


def make_filtering_merging_stats_plots(merging_stats):
    """Generate plotting dicts for merging statistics."""
    n_datasets = len(merging_stats)
    colors = [
        (color_list * int(math.ceil(n_datasets / len(color_list))))[i]
        for i in range(n_datasets)
    ]
    colors[-1] = "k"

    legends = ["initial scale"]
    if len(merging_stats) > 2:
        legends += [
            ordinal(i) + " rescale" for i in range(1, len(merging_stats) - 1)
        ] + ["final rescale"]
    elif len(merging_stats) == 2:
        legends += ["final rescale"]

    d = OrderedDict()
    overall_ccs = [m["overall"]["cc_one_half"] for m in merging_stats]
    overall_rpim = [m["overall"]["r_pim"] for m in merging_stats]
    overall_ioversigma = [m["overall"]["i_over_sigma_mean"] for m in merging_stats]
    overall_completeness = [m["overall"]["completeness"] for m in merging_stats]
    # first make overall plots
    d.update(
        {
            "cc_one_half_vs_cycle": {
                "data": [
                    {
                        "y": overall_ccs,
                        "x": range(1, n_datasets + 1),
                        "type": "scatter",
                        "mode": "lines",
                    }
                ],
                "layout": {
                    "title": "CC-half vs cycle",
                    "xaxis": {"title": "Cycle number"},
                    "yaxis": {"title": "CC-half"},
                },
            }
        }
    )
    d.update(
        {
            "r_pim_vs_cycle": {
                "data": [
                    {
                        "y": overall_rpim,
                        "x": range(1, n_datasets + 1),
                        "type": "scatter",
                        "mode": "lines",
                    }
                ],
                "layout": {
                    "title": "R-pim vs cycle",
                    "xaxis": {"title": "Cycle number"},
                    "yaxis": {"title": "R-pim"},
                },
            }
        }
    )
    d.update(
        {
            "i_over_sigma_vs_cycle": {
                "data": [
                    {
                        "y": overall_ioversigma,
                        "x": range(1, n_datasets + 1),
                        "type": "scatter",
                        "mode": "lines",
                    }
                ],
                "layout": {
                    "title": "<I/sigma> vs cycle",
                    "xaxis": {"title": "Cycle number"},
                    "yaxis": {"title": "<I/sigma>"},
                },
            }
        }
    )
    d.update(
        {
            "completeness_vs_cycle": {
                "data": [
                    {
                        "y": overall_completeness,
                        "x": range(1, n_datasets + 1),
                        "type": "scatter",
                        "mode": "lines",
                    }
                ],
                "layout": {
                    "title": "Completeness vs cycle",
                    "xaxis": {"title": "Cycle number"},
                    "yaxis": {"title": "Completeness"},
                },
            }
        }
    )
    cc_one_half_bins = merging_stats[0]["ccs"]
    r_pim_bins = merging_stats[0]["rpim"]
    r_merge_bins = merging_stats[0]["rmerge"]
    resolution = [1.0 / x ** 2 for x in merging_stats[0]["d_min"]]
    vals, txt = ResolutionPlotterMixin._d_star_sq_to_d_ticks(resolution, 5)
    d.update(
        {
            "cc_one_half_filter": {
                "data": [
                    {
                        "x": resolution,  # d_star_sq
                        "y": cc_one_half_bins,
                        "type": "scatter",
                        "name": legends[0],
                        "mode": "lines",
                        "line": {"color": colors[0]},
                    }
                ],
                "layout": {
                    "title": "CC-half vs resolution",
                    "xaxis": {
                        "title": u"Resolution (Å)",
                        "tickvals": vals,
                        "ticktext": txt,
                    },
                    "yaxis": {"title": "CC-half", "range": [0, 1]},
                },
            }
        }
    )
    d.update(
        {
            "r_pim_filter": {
                "data": [
                    {
                        "x": resolution,  # d_star_sq
                        "y": r_pim_bins,
                        "type": "scatter",
                        "name": legends[0],
                        "mode": "lines",
                        "line": {"color": colors[0]},
                    }
                ],
                "layout": {
                    "title": "R-pim vs resolution",
                    "xaxis": {
                        "title": u"Resolution (Å)",
                        "tickvals": vals,
                        "ticktext": txt,
                    },
                    "yaxis": {
                        "title": "R-pim",
                        "range": [0, min(1.5, max(r_pim_bins))],
                    },
                },
            }
        }
    )
    d.update(
        {
            "r_merge_filter": {
                "data": [
                    {
                        "x": resolution,  # d_star_sq
                        "y": r_merge_bins,
                        "type": "scatter",
                        "name": legends[0],
                        "mode": "lines",
                        "line": {"color": colors[0]},
                    }
                ],
                "layout": {
                    "title": "R-merge vs resolution",
                    "xaxis": {
                        "title": u"Resolution (Å)",
                        "tickvals": vals,
                        "ticktext": txt,
                    },
                    "yaxis": {
                        "title": "R-merge",
                        "range": [0, min(1.5, max(r_merge_bins))],
                    },
                },
            }
        }
    )
    for c, stats in enumerate(merging_stats[1:]):
        cc_one_half_bins = stats["ccs"]
        r_pim_bins = stats["rpim"]
        r_merge_bins = stats["rmerge"]
        resolution = [1.0 / x ** 2 for x in stats["d_min"]]
        d["cc_one_half_filter"]["data"].append(
            {
                "x": resolution,  # d_star_sq
                "y": cc_one_half_bins,
                "type": "scatter",
                "name": legends[c + 1],
                "line": {"color": colors[c + 1]},
            }
        )
        d["r_pim_filter"]["data"].append(
            {
                "x": resolution,  # d_star_sq
                "y": r_pim_bins,
                "type": "scatter",
                "name": legends[c + 1],
                "line": {"color": colors[c + 1]},
            }
        )
        d["r_merge_filter"]["data"].append(
            {
                "x": resolution,  # d_star_sq
                "y": r_merge_bins,
                "type": "scatter",
                "name": legends[c + 1],
                "line": {"color": colors[c + 1]},
            }
        )
    return d


def make_histogram_plots(cycle_results):
    """Make the histogram plots."""
    delta_cc_half_lists = [res["delta_cc_half_values"] for res in cycle_results]
    if not delta_cc_half_lists:
        return {}

    n = len(delta_cc_half_lists)

    d = OrderedDict()
    overall_mean_ccs = [res["mean_cc_half"] for res in cycle_results]
    d.update(
        {
            "mean_cc_one_half_vs_cycle": {
                "data": [
                    {
                        "y": overall_mean_ccs,
                        "x": range(1, n + 1),
                        "type": "scatter",
                        "mode": "lines",
                    }
                ],
                "layout": {
                    "title": "Resolution-averaged CC-half (sigma-tau) vs cycle",
                    "xaxis": {"title": "Cycle number"},
                    "yaxis": {"title": "Resolution-averaged CC-half (sigma-tau)"},
                },
            }
        }
    )
    colors = [(color_list * int(math.ceil(n / len(color_list))))[i] for i in range(n)]
    legends = [ordinal(i) + " Delta CC-half analysis" for i in range(1, n + 1)]
    if "image_ranges_removed" in cycle_results[0]:
        n_rej = [len(res["image_ranges_removed"]) for res in cycle_results]
    else:
        n_rej = [len(res["removed_datasets"]) for res in cycle_results]

    def _color_bar_charts(counts, index):
        n = 0
        bar_colors = []
        n_rej_this = n_rej[index]
        for count in counts:
            if n >= n_rej_this:
                bar_colors.append(colors[index])
            else:
                bar_colors.append("k")
                n += count
        return bar_colors

    def _add_new_histogram(d, hist, index):
        d.update(
            {
                "scale_filter_histograms_%s"
                % index: {
                    "data": [
                        {
                            "x": list(hist.slot_centers()),
                            "y": list(hist.slots()),
                            "type": "bar",
                            "name": legends[index],
                            "marker": {"color": _color_bar_charts(hist.slots(), index)},
                        }
                    ],
                    "layout": {
                        "title": "%s" % legends[index],
                        "xaxis": {"title": "Delta CC-half"},
                        "yaxis": {
                            "title": "Number of datasets/groups",
                            "range": [0, min(max(hist.slots()), 50)],
                        },
                        "bargap": 0,
                    },
                }
            }
        )

    for c, deltas in enumerate(delta_cc_half_lists):
        hist = flex.histogram(
            flex.double(deltas) * 100, min(deltas) * 100, max(deltas) * 100, n_slots=40
        )
        _add_new_histogram(d, hist, c)
    return d


def make_reduction_plots(initial_expids_and_image_ranges, expids_and_image_ranges):
    """Make a chart showing excluded image ranges."""
    x = list(range(len(initial_expids_and_image_ranges)))
    initial_n_images = []
    initial_expids = []
    for expid_and_img in initial_expids_and_image_ranges:
        initial_n_images.append(expid_and_img[1][1] - expid_and_img[1][0] + 1)
        initial_expids.append(expid_and_img[0])
    final_n_images = [0] * len(x)
    for expid, valid in expids_and_image_ranges:
        loc = initial_expids.index(expid)
        for t in valid:
            final_n_images[loc] = t[1] - t[0] + 1
    n_removed = [i - j for (i, j) in zip(initial_n_images, final_n_images)]
    d = {
        "reduction_plot": {
            "data": [
                {
                    "x": x,
                    "y": final_n_images,
                    "type": "bar",
                    "name": "final image ranges",
                },
                {"x": x, "y": n_removed, "type": "bar", "name": "removed image ranges"},
            ],
            "layout": {
                "title": "Image range plot",
                "xaxis": {"title": "Experiment number"},
                "yaxis": {"title": "Image number"},
                "bargap": 0,
                "barmode": "stack",
            },
        }
    }
    return d
