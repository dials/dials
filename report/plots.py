# -*- coding: utf-8 -*-
"""
This module defines a number of general plots, which may be relevant to
for reports of several programs.
"""
from __future__ import absolute_import, division, print_function
from collections import OrderedDict
import numpy as np
from cctbx import uctbx
from scitbx.array_family import flex
from scitbx.math import distributions


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


def i_over_sig_i_vs_i_plot(intensities, sigmas):
    """Plot unscaled I / sigma_adjusted vs unscaled I."""
    sel = (intensities > 0) & (sigmas > 0)
    intensities = intensities.select(sel)
    sigmas = sigmas.select(sel)
    x = flex.log10(intensities)
    y = intensities / sigmas

    H, xedges, yedges = np.histogram2d(
        x.as_numpy_array(), y.as_numpy_array(), bins=(200, 200)
    )
    nonzeros = np.nonzero(H)
    z = np.empty(H.shape)
    z[:] = np.NAN
    z[nonzeros] = H[nonzeros]

    y2 = intensities / (sigmas ** 2)
    H2, x2edges, y2edges = np.histogram2d(
        x.as_numpy_array(), y2.as_numpy_array(), bins=(200, 200)
    )
    nonzeros = np.nonzero(H2)
    z2 = np.empty(H2.shape)
    z2[:] = np.NAN
    z2[nonzeros] = H2[nonzeros]

    return {
        "i_over_sig_i_vs_i": {
            "data": [
                {
                    "x": xedges.tolist(),
                    "y": yedges.tolist(),
                    "z": z.transpose().tolist(),
                    "type": "heatmap",
                    "name": "Isigma distribution",
                    "colorbar": {
                        "title": "Number of reflections",
                        "titleside": "right",
                    },
                    "colorscale": "Jet",
                }
            ],
            "layout": {
                "title": "I/sig(I) vs I",
                "xaxis": {"title": "log I"},
                "yaxis": {"title": "I/sig(I)"},
            },
            "help": """\
This plot shows the distribution of I/sigma as a function of I, which can
give indication of the errors within the dataset. The I/sigma asymptotic
limit can be seen at the plateau in the top-right of the plot, if the measured
data are strong enough.

[1] Diederichs, K. (2010). Acta Cryst. D, 66(6), 733-740.
https://doi.org/10.1107/S0907444910014836
""",
        },
        "i_over_sig_isq_vs_i": {
            "data": [
                {
                    "x": x2edges.tolist(),
                    "y": y2edges.tolist(),
                    "z": z2.transpose().tolist(),
                    "type": "heatmap",
                    "name": "Isigma2 distribution",
                    "colorbar": {
                        "title": "Number of reflections",
                        "titleside": "right",
                    },
                    "colorscale": "Jet",
                }
            ],
            "layout": {
                "title": "I/sig(I)^2 vs I",
                "xaxis": {"title": "log I"},
                "yaxis": {"title": "I/sig(I)^2"},
            },
        },
    }


class ResolutionPlotterMixin(object):

    """Define additional helper methods for plotting"""

    @staticmethod
    def _d_star_sq_to_d_ticks(d_star_sq, nticks):
        min_d_star_sq = min(d_star_sq)
        dstep = (max(d_star_sq) - min_d_star_sq) / nticks
        tickvals = list(min_d_star_sq + (i * dstep) for i in range(nticks))
        ticktext = ["%.2f" % (uctbx.d_star_sq_as_d(dsq)) for dsq in tickvals]
        return tickvals, ticktext


class IntensityStatisticsPlots(ResolutionPlotterMixin):

    """Generate plots for intensity-derived statistics."""

    def __init__(
        self,
        intensities,
        anomalous=False,
        n_resolution_bins=20,
        xtriage_analyses=None,
        run_xtriage_analysis=True,
    ):
        self.n_bins = n_resolution_bins
        self._xanalysis = xtriage_analyses
        if anomalous:
            intensities = intensities.as_anomalous_array()
        intensities.setup_binner(n_bins=self.n_bins)
        merged = intensities.merge_equivalents()
        self.binner = intensities.binner()
        self.merged_intensities = merged.array()
        self.multiplicities = merged.redundancies().complete_array(new_data_value=0)
        if not self._xanalysis and run_xtriage_analysis:
            # imports needed here or won't work, unsure why.
            from mmtbx.scaling.xtriage import xtriage_analyses
            from mmtbx.scaling.xtriage import master_params as xtriage_master_params

            xtriage_params = xtriage_master_params.fetch(sources=[]).extract()
            xtriage_params.scaling.input.xray_data.skip_sanity_checks = True
            xanalysis = xtriage_analyses(
                miller_obs=self.merged_intensities,
                unmerged_obs=intensities,
                text_out="silent",
                params=xtriage_params,
            )
            self._xanalysis = xanalysis

    def generate_resolution_dependent_plots(self):
        d = OrderedDict()
        d.update(self.second_moments_plot())
        d.update(self.wilson_plot())
        return d

    def generate_miscellanous_plots(self):
        d = OrderedDict()
        d.update(self.cumulative_intensity_distribution_plot())
        d.update(self.l_test_plot())
        d.update(self.multiplicity_histogram())
        return d

    def multiplicity_histogram(self):
        """Generate histogram data for acentric and centric multiplicities."""
        mult_acentric = self.multiplicities.select_acentric().data()
        mult_centric = self.multiplicities.select_centric().data()

        multiplicities_acentric = {}
        multiplicities_centric = {}

        for x in sorted(set(mult_acentric)):
            multiplicities_acentric[x] = mult_acentric.count(x)
        for x in sorted(set(mult_centric)):
            multiplicities_centric[x] = mult_centric.count(x)

        return {
            "multiplicities": {
                "data": [
                    {
                        "x": list(multiplicities_acentric.keys()),
                        "y": list(multiplicities_acentric.values()),
                        "type": "bar",
                        "name": "Acentric",
                        "opacity": 0.75,
                    },
                    {
                        "x": list(multiplicities_centric.keys()),
                        "y": list(multiplicities_centric.values()),
                        "type": "bar",
                        "name": "Centric",
                        "opacity": 0.75,
                    },
                ],
                "layout": {
                    "title": "Distribution of multiplicities",
                    "xaxis": {"title": "Multiplicity"},
                    "yaxis": {
                        "title": "Frequency",
                        #'rangemode': 'tozero'
                    },
                    "bargap": 0,
                    "barmode": "overlay",
                },
            }
        }

    def wilson_plot(self):
        if not self._xanalysis or not self._xanalysis.wilson_scaling:
            return {}
        wilson_scaling = self._xanalysis.wilson_scaling
        tickvals_wilson, ticktext_wilson = self._d_star_sq_to_d_ticks(
            wilson_scaling.d_star_sq, nticks=5
        )

        return {
            "wilson_intensity_plot": {
                "data": (
                    [
                        {
                            "x": list(wilson_scaling.d_star_sq),
                            "y": list(wilson_scaling.mean_I_obs_data),
                            "type": "scatter",
                            "name": "Observed",
                        },
                        {
                            "x": list(wilson_scaling.d_star_sq),
                            "y": list(wilson_scaling.mean_I_obs_theory),
                            "type": "scatter",
                            "name": "Expected",
                        },
                        {
                            "x": list(wilson_scaling.d_star_sq),
                            "y": list(wilson_scaling.mean_I_normalisation),
                            "type": "scatter",
                            "name": "Smoothed",
                        },
                    ]
                ),
                "layout": {
                    "title": "Wilson intensity plot",
                    "xaxis": {
                        "title": u"Resolution (Å)",
                        "tickvals": tickvals_wilson,
                        "ticktext": ticktext_wilson,
                    },
                    "yaxis": {"type": "log", "title": "Mean(I)", "rangemode": "tozero"},
                },
            }
        }

    def cumulative_intensity_distribution_plot(self):
        if not self._xanalysis or not self._xanalysis.twin_results:
            return {}
        nz_test = self._xanalysis.twin_results.nz_test
        return {
            "cumulative_intensity_distribution": {
                "data": [
                    {
                        "x": list(nz_test.z),
                        "y": list(nz_test.ac_obs),
                        "type": "scatter",
                        "name": "Acentric observed",
                        "mode": "lines",
                        "line": {"color": "rgb(31, 119, 180)"},
                    },
                    {
                        "x": list(nz_test.z),
                        "y": list(nz_test.c_obs),
                        "type": "scatter",
                        "name": "Centric observed",
                        "mode": "lines",
                        "line": {"color": "rgb(255, 127, 14)"},
                    },
                    {
                        "x": list(nz_test.z),
                        "y": list(nz_test.ac_untwinned),
                        "type": "scatter",
                        "name": "Acentric theory",
                        "mode": "lines",
                        "line": {"color": "rgb(31, 119, 180)", "dash": "dot"},
                        "opacity": 0.8,
                    },
                    {
                        "x": list(nz_test.z),
                        "y": list(nz_test.c_untwinned),
                        "type": "scatter",
                        "name": "Centric theory",
                        "mode": "lines",
                        "line": {"color": "rgb(255, 127, 14)", "dash": "dot"},
                        "opacity": 0.8,
                    },
                ],
                "layout": {
                    "title": "Cumulative intensity distribution",
                    "xaxis": {"title": "z", "range": (0, 1)},
                    "yaxis": {"title": "P(Z <= Z)", "range": (0, 1)},
                },
            }
        }

    def l_test_plot(self):
        if not self._xanalysis or not self._xanalysis.twin_results:
            return {}
        l_test = self._xanalysis.twin_results.l_test
        return {
            "l_test": {
                "data": [
                    {
                        "x": list(l_test.l_values),
                        "y": list(l_test.l_cumul_untwinned),
                        "type": "scatter",
                        "name": "Untwinned",
                        "mode": "lines",
                        "line": {"color": "rgb(31, 119, 180)", "dash": "dashdot"},
                    },
                    {
                        "x": list(l_test.l_values),
                        "y": list(l_test.l_cumul_perfect_twin),
                        "type": "scatter",
                        "name": "Perfect twin",
                        "mode": "lines",
                        "line": {"color": "rgb(31, 119, 180)", "dash": "dot"},
                        "opacity": 0.8,
                    },
                    {
                        "x": list(l_test.l_values),
                        "y": list(l_test.l_cumul),
                        "type": "scatter",
                        "name": "Observed",
                        "mode": "lines",
                        "line": {"color": "rgb(255, 127, 14)"},
                    },
                ],
                "layout": {
                    "title": "L test (Padilla and Yeates)",
                    "xaxis": {"title": "|l|", "range": (0, 1)},
                    "yaxis": {"title": "P(L >= l)", "range": (0, 1)},
                },
            }
        }

    def second_moments_plot(self):

        acentric = self.merged_intensities.select_acentric()
        centric = self.merged_intensities.select_centric()
        if acentric.size():
            acentric.setup_binner(n_bins=self.n_bins)
            second_moments_acentric = acentric.second_moment_of_intensities(
                use_binning=True
            )
        else:
            second_moments_acentric = None
        if centric.size():
            centric.setup_binner(n_bins=self.n_bins)
            second_moments_centric = centric.second_moment_of_intensities(
                use_binning=True
            )
        else:
            second_moments_centric = None

        second_moment_d_star_sq = []
        if acentric.size():
            second_moment_d_star_sq.extend(
                second_moments_acentric.binner.bin_centers(2)
            )
        if centric.size():
            second_moment_d_star_sq.extend(second_moments_centric.binner.bin_centers(2))
        tickvals_2nd_moment, ticktext_2nd_moment = self._d_star_sq_to_d_ticks(
            second_moment_d_star_sq, nticks=5
        )

        return {
            "second_moments": {
                "data": [
                    (
                        {
                            "x": list(
                                second_moments_acentric.binner.bin_centers(2)
                            ),  # d_star_sq
                            "y": second_moments_acentric.data[1:-1],
                            "type": "scatter",
                            "name": "<I^2> acentric",
                        }
                        if acentric.size()
                        else {}
                    ),
                    (
                        {
                            "x": list(
                                second_moments_centric.binner.bin_centers(2)
                            ),  # d_star_sq
                            "y": second_moments_centric.data[1:-1],
                            "type": "scatter",
                            "name": "<I^2> centric",
                        }
                        if centric.size()
                        else {}
                    ),
                ],
                "layout": {
                    "title": "Second moment of I",
                    "xaxis": {
                        "title": u"Resolution (Å)",
                        "tickvals": tickvals_2nd_moment,
                        "ticktext": ticktext_2nd_moment,
                    },
                    "yaxis": {"title": "<I^2>", "rangemode": "tozero"},
                },
            }
        }


class ResolutionPlotsAndStats(ResolutionPlotterMixin):

    """
    Use iotbx dataset statistics objects to make plots and tables for reports.

    This class allows the generation of plots of various properties as a
    function of resolution as well as a statistics table and summary table,
    using the data from two iotbx.dataset_statistics objects, with
    anomalous=False/True.
    """

    def __init__(
        self, dataset_statistics, anomalous_dataset_statistics, is_centric=False
    ):
        self.dataset_statistics = dataset_statistics
        self.anomalous_dataset_statistics = anomalous_dataset_statistics
        self.d_star_sq_bins = [
            (1 / bin_stats.d_min ** 2) for bin_stats in self.dataset_statistics.bins
        ]
        self.d_star_sq_tickvals, self.d_star_sq_ticktext = self._d_star_sq_to_d_ticks(
            self.d_star_sq_bins, nticks=5
        )
        self.is_centric = is_centric

    def make_all_plots(self):
        """Make a dictionary containing all available resolution-dependent plots."""
        d = OrderedDict()
        d.update(self.cc_one_half_plot())
        d.update(self.i_over_sig_i_plot())
        d.update(self.completeness_plot())
        d.update(self.multiplicity_vs_resolution_plot())
        return d

    def cc_one_half_plot(self, method=None):
        """Make a plot of cc half against resolution."""

        if method == "sigma_tau":
            cc_one_half_bins = [
                bin_stats.cc_one_half_sigma_tau
                if bin_stats.cc_one_half_sigma_tau
                else 0.0
                for bin_stats in self.dataset_statistics.bins
            ]
            cc_one_half_critical_value_bins = [
                bin_stats.cc_one_half_sigma_tau_critical_value
                if bin_stats.cc_one_half_sigma_tau_critical_value
                else 0.0
                for bin_stats in self.dataset_statistics.bins
            ]
        else:
            cc_one_half_bins = [
                bin_stats.cc_one_half if bin_stats.cc_one_half else 0.0
                for bin_stats in self.dataset_statistics.bins
            ]
            cc_one_half_critical_value_bins = [
                bin_stats.cc_one_half_critical_value
                if bin_stats.cc_one_half_critical_value
                else 0.0
                for bin_stats in self.dataset_statistics.bins
            ]
        cc_anom_bins = [
            bin_stats.cc_anom if bin_stats.cc_anom else 0.0
            for bin_stats in self.dataset_statistics.bins
        ]
        cc_anom_critical_value_bins = [
            bin_stats.cc_anom_critical_value
            if bin_stats.cc_anom_critical_value
            else 0.0
            for bin_stats in self.dataset_statistics.bins
        ]

        return {
            "cc_one_half": {
                "data": [
                    {
                        "x": self.d_star_sq_bins,  # d_star_sq
                        "y": cc_one_half_bins,
                        "type": "scatter",
                        "name": "CC-half",
                        "mode": "lines",
                        "line": {"color": "rgb(31, 119, 180)"},
                    },
                    {
                        "x": self.d_star_sq_bins,  # d_star_sq
                        "y": cc_one_half_critical_value_bins,
                        "type": "scatter",
                        "name": "CC-half critical value (p=0.01)",
                        "line": {"color": "rgb(31, 119, 180)", "dash": "dot"},
                    },
                    (
                        {
                            "x": self.d_star_sq_bins,  # d_star_sq
                            "y": cc_anom_bins,
                            "type": "scatter",
                            "name": "CC-anom",
                            "mode": "lines",
                            "line": {"color": "rgb(255, 127, 14)"},
                        }
                        if not self.is_centric
                        else {}
                    ),
                    (
                        {
                            "x": self.d_star_sq_bins,  # d_star_sq
                            "y": cc_anom_critical_value_bins,
                            "type": "scatter",
                            "name": "CC-anom critical value (p=0.01)",
                            "mode": "lines",
                            "line": {"color": "rgb(255, 127, 14)", "dash": "dot"},
                        }
                        if not self.is_centric
                        else {}
                    ),
                ],
                "layout": {
                    "title": "CC-half vs resolution",
                    "xaxis": {
                        "title": u"Resolution (Å)",
                        "tickvals": self.d_star_sq_tickvals,
                        "ticktext": self.d_star_sq_ticktext,
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

    def i_over_sig_i_plot(self):
        """Make a plot of <I/sigI> against resolution."""
        i_over_sig_i_bins = [
            bin_stats.i_over_sigma_mean for bin_stats in self.dataset_statistics.bins
        ]

        return {
            "i_over_sig_i": {
                "data": [
                    {
                        "x": self.d_star_sq_bins,  # d_star_sq
                        "y": i_over_sig_i_bins,
                        "type": "scatter",
                        "name": "I/sigI vs resolution",
                    }
                ],
                "layout": {
                    "title": "<I/sig(I)> vs resolution",
                    "xaxis": {
                        "title": u"Resolution (Å)",
                        "tickvals": self.d_star_sq_tickvals,
                        "ticktext": self.d_star_sq_ticktext,
                    },
                    "yaxis": {"title": "<I/sig(I)>", "rangemode": "tozero"},
                },
            }
        }

    def completeness_plot(self):
        """Make a plot of completeness against resolution."""
        completeness_bins = [
            bin_stats.completeness for bin_stats in self.dataset_statistics.bins
        ]
        anom_completeness_bins = [
            bin_stats.anom_completeness
            for bin_stats in self.anomalous_dataset_statistics.bins
        ]

        return {
            "completeness": {
                "data": [
                    {
                        "x": self.d_star_sq_bins,
                        "y": completeness_bins,
                        "type": "scatter",
                        "name": "Completeness",
                    },
                    (
                        {
                            "x": self.d_star_sq_bins,
                            "y": anom_completeness_bins,
                            "type": "scatter",
                            "name": "Anomalous completeness",
                        }
                        if not self.is_centric
                        else {}
                    ),
                ],
                "layout": {
                    "title": "Completeness vs resolution",
                    "xaxis": {
                        "title": u"Resolution (Å)",
                        "tickvals": self.d_star_sq_tickvals,
                        "ticktext": self.d_star_sq_ticktext,
                    },
                    "yaxis": {"title": "Completeness", "range": (0, 1)},
                },
            }
        }

    def multiplicity_vs_resolution_plot(self):
        """Make a plot of multiplicity against resolution."""
        multiplicity_bins = [
            bin_stats.mean_redundancy for bin_stats in self.dataset_statistics.bins
        ]
        anom_multiplicity_bins = [
            bin_stats.mean_redundancy
            for bin_stats in self.anomalous_dataset_statistics.bins
        ]

        return {
            "multiplicity_vs_resolution": {
                "data": [
                    {
                        "x": self.d_star_sq_bins,
                        "y": multiplicity_bins,
                        "type": "scatter",
                        "name": "Multiplicity",
                    },
                    (
                        {
                            "x": self.d_star_sq_bins,
                            "y": anom_multiplicity_bins,
                            "type": "scatter",
                            "name": "Anomalous multiplicity",
                        }
                        if not self.is_centric
                        else {}
                    ),
                ],
                "layout": {
                    "title": "Multiplicity vs resolution",
                    "xaxis": {
                        "title": u"Resolution (Å)",
                        "tickvals": self.d_star_sq_tickvals,
                        "ticktext": self.d_star_sq_ticktext,
                    },
                    "yaxis": {"title": "Multiplicity", "rangemode": "tozero"},
                },
            }
        }

    def merging_statistics_table(self, cc_half_method=None):

        headers = [
            u"Resolution (Å)",
            "N(obs)",
            "N(unique)",
            "Multiplicity",
            "Completeness",
            "Mean(I)",
            "Mean(I/sigma)",
            "Rmerge",
            "Rmeas",
            "Rpim",
            "CC1/2",
        ]
        if not self.is_centric:
            headers.append("CCano")
        rows = []

        def safe_format(format_str, item):
            return format_str % item if item is not None else ""

        for bin_stats in self.dataset_statistics.bins:
            row = [
                "%.2f - %.2f" % (bin_stats.d_max, bin_stats.d_min),
                bin_stats.n_obs,
                bin_stats.n_uniq,
                "%.2f" % bin_stats.mean_redundancy,
                "%.2f" % (100 * bin_stats.completeness),
                "%.1f" % bin_stats.i_mean,
                "%.1f" % bin_stats.i_over_sigma_mean,
                safe_format("%.3f", bin_stats.r_merge),
                safe_format("%.3f", bin_stats.r_meas),
                safe_format("%.3f", bin_stats.r_pim),
            ]
            if cc_half_method == "sigma_tau":
                row.append(
                    "%.3f%s"
                    % (
                        bin_stats.cc_one_half_sigma_tau,
                        "*" if bin_stats.cc_one_half_sigma_tau_significance else "",
                    )
                )
            else:
                row.append(
                    "%.3f%s"
                    % (
                        bin_stats.cc_one_half,
                        "*" if bin_stats.cc_one_half_significance else "",
                    )
                )

            if not self.is_centric:
                row.append(
                    "%.3f%s"
                    % (bin_stats.cc_anom, "*" if bin_stats.cc_anom_significance else "")
                )
            rows.append(row)

        merging_stats_table = [headers]
        merging_stats_table.extend(rows)

        return merging_stats_table

    def overall_statistics_table(self, cc_half_method=None):

        headers = ["", "Overall", "Low resolution", "High resolution"]

        stats = (
            self.dataset_statistics.overall,
            self.dataset_statistics.bins[0],
            self.dataset_statistics.bins[-1],
        )

        rows = [
            [u"Resolution (Å)"] + ["%.2f - %.2f" % (s.d_max, s.d_min) for s in stats],
            ["Observations"] + ["%i" % s.n_obs for s in stats],
            ["Unique reflections"] + ["%i" % s.n_uniq for s in stats],
            ["Multiplicity"] + ["%.1f" % s.mean_redundancy for s in stats],
            ["Completeness"] + ["%.2f%%" % (s.completeness * 100) for s in stats],
            # ['Mean intensity'] + ['%.1f' %s.i_mean for s in stats],
            ["Mean I/sigma(I)"] + ["%.1f" % s.i_over_sigma_mean for s in stats],
            ["Rmerge"] + ["%.3f" % s.r_merge for s in stats],
            ["Rmeas"] + ["%.3f" % s.r_meas for s in stats],
            ["Rpim"] + ["%.3f" % s.r_pim for s in stats],
        ]

        if cc_half_method == "sigma_tau":
            rows.append(["CC1/2"] + ["%.3f" % s.cc_one_half_sigma_tau for s in stats])
        else:
            rows.append(["CC1/2"] + ["%.3f" % s.cc_one_half for s in stats])
        rows = [[u"<strong>%s</strong>" % r[0]] + r[1:] for r in rows]

        overall_stats_table = [headers]
        overall_stats_table.extend(rows)

        return overall_stats_table

    def statistics_tables(self):
        """Generate the overall and by-resolution tables."""
        return (self.overall_statistics_table(), self.merging_statistics_table())


class AnomalousPlotter(ResolutionPlotterMixin):
    def __init__(self, anomalous_array, strong_cutoff=0.0, n_bins=20):
        self.intensities_anom = anomalous_array.map_to_asu()
        self.merged = self.intensities_anom.merge_equivalents(
            use_internal_variance=False
        ).array()
        self.n_bins = n_bins
        self.strong_cutoff = strong_cutoff
        if strong_cutoff > 0.0:
            self.low_res_intensities_anom = self.intensities_anom.resolution_filter(
                d_min=strong_cutoff
            )
            self.strong_merged = self.low_res_intensities_anom.merge_equivalents(
                use_internal_variance=False
            ).array()

    def make_plots(self):
        d = OrderedDict()
        if self.strong_cutoff > 0.0:
            d.update(self.del_anom_normal_plot(self.strong_merged, self.strong_cutoff))
            d.update(
                self.del_anom_scatter_plot(
                    self.low_res_intensities_anom, self.strong_cutoff
                )
            )
        else:
            d.update(self.del_anom_normal_plot(self.merged))
        d.update(self.del_anom_correlation_ratio(self.intensities_anom))
        return d

    def del_anom_correlation_ratio(self, unmerged_intensities):

        acentric = unmerged_intensities.select_acentric()
        centric = unmerged_intensities.select_centric()
        correl_ratios_acentric, correl_ratios_centric = ([], [])

        def calc_correl_ratios(data):
            correl_ratios = []
            data.setup_binner(n_bins=self.n_bins)
            for i_bin in data.binner().range_used():
                sel = data.binner().selection(i_bin)
                data_sel = data.select(sel)
                if data_sel.size() > 0:
                    arr1, arr2 = data_sel.half_dataset_anomalous_correlation(
                        return_split_datasets=1
                    )
                    dano1 = arr1.anomalous_differences().data()
                    dano2 = arr2.anomalous_differences().data()
                    if dano1.size() > 0:
                        rmsd_11 = (
                            flex.sum((dano1 - dano2) ** 2) / (2.0 * dano1.size())
                        ) ** 0.5
                        rmsd_1min1 = (
                            flex.sum((dano1 + dano2) ** 2) / (2.0 * dano1.size())
                        ) ** 0.5
                        correl_ratios.append(rmsd_1min1 / rmsd_11)
                    else:
                        correl_ratios.append(0.0)
                else:
                    correl_ratios.append(0.0)
            return correl_ratios

        if acentric.size() > 0:
            correl_ratios_acentric = calc_correl_ratios(acentric)
            if all(list(flex.double(correl_ratios_acentric) == 0.0)):
                correl_ratios_acentric = []
            else:
                d_star_sq_acentric = acentric.binner().bin_centers(2)
                actickvals, acticktext = self._d_star_sq_to_d_ticks(
                    d_star_sq_acentric, nticks=5
                )
        if centric.size() > 0:
            correl_ratios_centric = calc_correl_ratios(centric)
            if all(list(flex.double(correl_ratios_centric) == 0.0)):
                correl_ratios_centric = []
            else:
                d_star_sq_centric = centric.binner().bin_centers(2)
                ctickvals, cticktext = self._d_star_sq_to_d_ticks(
                    d_star_sq_acentric, nticks=5
                )

        if not (correl_ratios_acentric or correl_ratios_centric):
            return {}
        if correl_ratios_acentric:
            tickvals = actickvals
            ticktext = acticktext
        else:
            tickvals = ctickvals
            ticktext = cticktext
        return {
            "anom_correl_plot": {
                "data": [
                    (
                        {
                            "x": list(d_star_sq_acentric),
                            "y": correl_ratios_acentric,
                            "type": "lines",
                            "name": "Anomalous correlation ratio (acentric)",
                        }
                        if correl_ratios_acentric
                        else {}
                    ),
                    (
                        {
                            "x": list(d_star_sq_centric),
                            "y": correl_ratios_centric,
                            "type": "lines",
                            "name": "Anomalous correlation ratio (centric)",
                        }
                        if correl_ratios_centric
                        else {}
                    ),
                ],
                "layout": {
                    "title": "Anomalous R.M.S. correlation ratio (acentric reflections)",
                    "xaxis": {
                        "title": u"Resolution (Å)",
                        "tickvals": tickvals,
                        "ticktext": ticktext,
                    },
                    "yaxis": {"anchor": "x", "title": "rms correlation ratio"},
                },
                "help": """\
This plot shows the significance of the anomalous signal, as shown in the
anomalous scatter plot, by calculating the ratio of the width of the signal along
the diagonal (a measure of the anomalous signal) over the width of the signal
perpendicular to the diagonal (a measure of the error).

[1] P. Evans, Acta Cryst. (2006). D62, 72-82
https://doi.org/10.1107/S0907444905036693
""",
            }
        }

    def del_anom_scatter_plot(self, unmerged_intensities, strong_cutoff=0.0):
        """Make a scatter plot of the anomalous differences of half sets."""

        acentric = unmerged_intensities.select_acentric()
        if acentric.size() == 0:
            return {}
        arr1, arr2 = acentric.half_dataset_anomalous_correlation(
            return_split_datasets=1
        )
        dano1 = arr1.anomalous_differences()
        dano2 = arr2.anomalous_differences()
        assert dano1.indices().all_eq(dano2.indices())
        if dano1.size() == 0:
            return {}
        max_val = max(flex.max(dano1.data()), flex.max(dano2.data()))
        min_val = min(flex.min(dano1.data()), flex.min(dano2.data()))

        title = "Correlation of half-set differences"
        plotname = "anom_scatter_plot"
        if strong_cutoff > 0.0:
            title += " (d > %.2f)" % strong_cutoff
            plotname += "_lowres"
        else:
            title += " (all data)"
        return {
            plotname: {
                "data": [
                    {
                        "x": list(dano1.data()),
                        "y": list(dano2.data()),
                        "type": "scatter",
                        "mode": "markers",
                        "size": 1,
                        "name": "half-set anomalous differences (acentrics)",
                    },
                    {
                        "x": [min_val - 1, max_val + 1],
                        "y": [min_val - 1, max_val + 1],
                        "type": "scatter",
                        "mode": "lines",
                        "name": "D1 = D2",
                        "color": "rgb(0,0,0)",
                    },
                ],
                "layout": {
                    "title": title,
                    "xaxis": {"anchor": "y", "title": "Delta I1"},
                    "yaxis": {"anchor": "x", "title": "Delta I2"},
                },
                "help": """\
This plot shows the correlation of the anomalous differences for the data divided
into two half sets. For each reflection, the I+ and I- observations are divided
into two sets, and two differences are calculated; Delta I1 = I+(1) - I-(1),
Delta I2 = I+(2) - I-(2). Perfect data would therefore have all points along
the diagonal, in reality an elliptical distribution is seen in the presence of
anomalous signal, or a spherical distribution for data with no anomalous signal.

[1] P. Evans, Acta Cryst. (2006). D62, 72-82
https://doi.org/10.1107/S0907444905036693
""",
            }
        }

    @staticmethod
    def del_anom_normal_plot(intensities, strong_cutoff=0.0):
        """Make a normal probability plot of the normalised anomalous differences."""
        diff_array = intensities.anomalous_differences()
        if not diff_array.data().size():
            return {}
        delta = diff_array.data() / diff_array.sigmas()

        norm = distributions.normal_distribution()

        n = len(delta)
        if n <= 10:
            a = 3 / 8
        else:
            a = 0.5

        y = flex.sorted(delta)
        x = [norm.quantile((i + 1 - a) / (n + 1 - (2 * a))) for i in range(n)]

        H, xedges, yedges = np.histogram2d(
            np.array(x), y.as_numpy_array(), bins=(200, 200)
        )
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

        title = "Normal probability plot of anomalous differences"
        plotname = "normal_distribution_plot"
        if strong_cutoff > 0.0:
            title += " (d > %.2f)" % strong_cutoff
            plotname += "_lowres"
        else:
            title += " (all data)"
            plotname += "_highres"
        return {
            plotname: {
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
                    "title": title,
                    "xaxis": {
                        "anchor": "y",
                        "title": "expected delta",
                        "range": [-4, 4],
                    },
                    "yaxis": {
                        "anchor": "x",
                        "title": "observed delta",
                        "range": [-5, 5],
                    },
                },
                "help": """\
    This plot shows the normalised anomalous differences, sorted in order and
    plotted against the expected order based on a normal distribution model.
    A true normal distribution of deviations would give the straight line indicated.

    [1] P. L. Howell and G. D. Smith, J. Appl. Cryst. (1992). 25, 81-86
    https://doi.org/10.1107/S0021889891010385
    [2] P. Evans, Acta Cryst. (2006). D62, 72-82
    https://doi.org/10.1107/S0907444905036693
    """,
            }
        }
