"""
Algorithms for analysis of resolution limits.
"""

from __future__ import annotations

import enum
import logging
import math
import typing

import iotbx.merging_statistics
import iotbx.mtz
import iotbx.phil
from cctbx import miller, uctbx
from cctbx.array_family import flex
from iotbx.reflection_file_utils import label_table
from scitbx.math import curve_fitting, five_number_summary

from dials.algorithms.scaling.scaling_library import determine_best_unit_cell
from dials.report import plots
from dials.util import Sorry, tabulate
from dials.util.batch_handling import (
    assign_batches_to_reflections,
    calculate_batch_offsets,
)
from dials.util.filter_reflections import filter_reflection_table
from dials.util.normalisation import quasi_normalisation

logger = logging.getLogger(__name__)


class metrics(enum.Enum):
    """Supported metrics for estimation of resolution limits."""

    CC_HALF = "cc_half"
    CC_REF = "cc_ref"
    ISIGMA = "unmerged_i_over_sigma_mean"
    MISIGMA = "i_over_sigma_mean"
    I_MEAN_OVER_SIGMA_MEAN = "i_mean_over_sigi_mean"
    RMERGE = "r_merge"
    COMPLETENESS = "completeness"


def polynomial_fit(x, y, degree=5):
    """
    Fit a polynomial to the values y(x) and return this fit

    x, y should be iterables containing floats of the same size. The order is the order
    of polynomial to use for this fit. This will be useful for e.g. I/sigma.
    """

    fit = curve_fitting.univariate_polynomial_fit(
        x, y, degree=degree, max_iterations=100
    )
    f = curve_fitting.univariate_polynomial(*fit.params)
    return f(x)


def tanh_fit(x, y, iqr_multiplier=None):
    """
    Fit a tanh function to the values y(x) and return this fit

    x, y should be iterables containing floats of the same size. This is used for
    fitting a curve to CC½.
    """

    tf = curve_fitting.tanh_fit(x, y)
    f = curve_fitting.tanh(*tf.params)

    if iqr_multiplier:
        assert iqr_multiplier > 0
        yc = f(x)
        dy = y - yc

        min_x, q1_x, med_x, q3_x, max_x = five_number_summary(dy)
        iqr_x = q3_x - q1_x
        cut_x = iqr_multiplier * iqr_x
        outliers = (dy > q3_x + cut_x) | (dy < q1_x - cut_x)
        if outliers.count(True) > 0:
            xo = x.select(~outliers)
            yo = y.select(~outliers)
            tf = curve_fitting.tanh_fit(xo, yo)
            f = curve_fitting.tanh(*tf.params)

    return f(x)


def log_fit(x, y, degree=5):
    """Fit the values log(y(x)) then return exp() to this fit.

    x, y should be iterables containing floats of the same size. The order is the order
    of polynomial to use for this fit. This will be useful for e.g. I/sigma."""

    fit = curve_fitting.univariate_polynomial_fit(
        x, flex.log(y), degree=degree, max_iterations=100
    )
    f = curve_fitting.univariate_polynomial(*fit.params)
    return flex.exp(f(x))


def log_inv_fit(x, y, degree=5):
    """Fit the values log(1 / y(x)) then return the inverse of this fit.

    x, y should be iterables, the order of the polynomial for the transformed
    fit needs to be specified. This will be useful for e.g. Rmerge."""

    fit = curve_fitting.univariate_polynomial_fit(
        x, flex.log(1 / y), degree=degree, max_iterations=100
    )
    f = curve_fitting.univariate_polynomial(*fit.params)
    return 1 / flex.exp(f(x))


def resolution_fit_from_merging_stats(merging_stats, metric, model, limit, sel=None):
    """Estimate a resolution limit based on the input `metric`

    The function defined by `model` will be fit to the selected `metric` which has been
    pre-calculated by the `merging_stats` object. The estimated resolution limit is
    chosen as the `d_star_sq` value at which the fitted function equals `limit`.

    Args:
        merging_stats (iotbx.merging_statistics.dataset_statistics): Pre-calculated
            merging statistics object
        metric (str): The metric to use for estimating a resolution limit. Must be a
            metric calculated by `iotbx.merging_statistics.merging_stats` and
            available as an attribute on the `bins` attribute of the input
            `merging_stats` object.
        model: The function to fit to the selected metric. Must be callable, taking as
            input x (d_star_sq) and y (the metric to be fitted) values, returning the
            fitted y(x) values.
        limit (float): The resolution limit criterion.
        sel (scitbx.array_family.flex.bool): An optional selection to apply to the
            `merging_stats` bins.

    Returns: The estimated resolution limit in units of Å^-1
    """

    y_obs = flex.double(getattr(b, metric) for b in merging_stats.bins).reversed()
    d_star_sq = flex.double(
        uctbx.d_as_d_star_sq(b.d_min) for b in merging_stats.bins
    ).reversed()
    return resolution_fit(d_star_sq, y_obs, model, limit, sel=sel)


def resolution_fit(d_star_sq, y_obs, model, limit, sel=None):
    """Estimate a resolution limit based on the input merging statistics

    The function defined by `model` will be fit to the input `d_star_sq` and `y_obs`.
    The estimated resolution limit is chosen as the `d_star_sq` value at which the
    fitted function equals `limit`.

    Args:
        d_star_sq (scitbx.array_family.flex.double): The high resolution limits of the
            resolution bins in units 1/d*2
        y_obs (scitbx.array_family.flex.double): The statistic against which to fit the
            function `model`
        model: The function to fit against `y_obs`. Must be callable, taking as input x
            (d_star_sq) and y (the metric to be fitted) values, returning the fitted
            y(x) values.
        limit (float): The resolution limit criterion.
        sel (scitbx.array_family.flex.bool): An optional selection to apply to the
            `d_star_sq` and `y_obs` values.

    Returns: The estimated resolution limit in units of Å^-1

    Raises:
        RuntimeError: Raised if no `y_obs` values remain after application of the
        selection `sel`
    """
    if not sel:
        sel = flex.bool(len(d_star_sq), True)
    sel &= y_obs > 0
    y_obs = y_obs.select(sel)
    d_star_sq = d_star_sq.select(sel)

    if not len(y_obs):
        raise RuntimeError("No reflections left for fitting")
    y_fit = model(d_star_sq, y_obs, 6)
    logger.debug(
        tabulate(
            [("d*2", "d", "obs", "fit")]
            + [
                (ds2, uctbx.d_star_sq_as_d(ds2), yo, yf)
                for ds2, yo, yf in zip(d_star_sq, y_obs, y_fit)
            ],
            headers="firstrow",
        )
    )

    if flex.min(y_obs) > limit:
        d_min = 1.0 / math.sqrt(flex.max(d_star_sq))
    else:
        try:
            d_min = 1.0 / math.sqrt(interpolate_value(d_star_sq, y_fit, limit))
        except RuntimeError as e:
            logger.debug(f"Error interpolating value: {e}")
            d_min = None

    return ResolutionResult(d_star_sq, y_obs, y_fit, d_min)


def _get_cc_half_significance(merging_stats, cc_half_method):
    """Get the CC½ significance values from the input merging_stats object"""
    if (
        cc_half_method == "sigma_tau"
        and merging_stats.overall.cc_one_half_sigma_tau_significance is not None
    ):
        return flex.bool(
            b.cc_one_half_sigma_tau_significance for b in merging_stats.bins
        ).reversed()
    elif merging_stats.overall.cc_one_half_significance is not None:
        return flex.bool(
            b.cc_one_half_significance for b in merging_stats.bins
        ).reversed()


def _get_cc_half_critical_values(merging_stats, cc_half_method):
    """Get the CC½ critical values from the input merging_stats object"""
    if (
        cc_half_method == "sigma_tau"
        and merging_stats.overall.cc_one_half_sigma_tau_critical_value is not None
    ):
        return flex.double(
            b.cc_one_half_sigma_tau_critical_value for b in merging_stats.bins
        ).reversed()
    elif merging_stats.overall.cc_one_half_critical_value is not None:
        critical = [
            b.cc_one_half_critical_value
            if b.cc_one_half_critical_value is not None
            else 0.0
            for b in merging_stats.bins
        ]
        return flex.double(critical).reversed()


def resolution_cc_half(
    merging_stats, limit, cc_half_method="half_dataset", model=tanh_fit
):
    """Estimate a resolution limit based on CC½

    The function defined by `model` will be fit to the CC½ values that have been
    pre-calculated by the `merging_stats` object. The estimated resolution limit is
    chosen as the `d_star_sq` value at which the fitted function equals `limit`.

    Args:
        merging_stats (iotbx.merging_statistics.dataset_statistics): Pre-calculated
            merging statistics object
        cc_half_method (str): The method for calculating CC½. Either "half_dataset" or
            "sigma_tau" (See Assmann et al., J. Appl. Cryst. (2016). 49, 1021–1028).
        model: The function to fit to the selected metric. Must be callable, taking as
            input x (d_star_sq) and y (the metric to be fitted) values, returning the
            fitted y(x) values. Default is `tanh_fit`.
        limit (float): The resolution limit criterion.

    Returns: The estimated resolution limit in units of Å^-1
    """
    sel = _get_cc_half_significance(merging_stats, cc_half_method)
    metric = "cc_one_half_sigma_tau" if cc_half_method == "sigma_tau" else "cc_one_half"
    result = resolution_fit_from_merging_stats(
        merging_stats, metric, model, limit, sel=sel
    )
    critical_values = _get_cc_half_critical_values(merging_stats, cc_half_method)
    if critical_values:
        result = result._replace(critical_values=critical_values.select(sel))
    return result


def interpolate_value(x, y, t):
    """Find the value of x: y(x) = t."""

    if t > max(y) or t < min(y):
        raise RuntimeError(f"t outside of [{min(y):f}, {max(y):f}]")

    for j in range(1, len(x)):
        x0 = x[j - 1]
        y0 = y[j - 1]

        x1 = x[j]
        y1 = y[j]

        if (y0 - t) * (y1 - t) < 0:
            return x0 + (t - y0) * (x1 - x0) / (y1 - y0)


def miller_array_from_mtz(unmerged_mtz, anomalous=False, labels=None):
    mtz_object = iotbx.mtz.object(file_name=unmerged_mtz)
    miller_arrays = mtz_object.as_miller_arrays(
        merge_equivalents=False, anomalous=anomalous
    )
    i_obs = None
    batches = None
    all_i_obs = []
    for array in miller_arrays:
        labels = array.info().label_string()
        if array.is_xray_intensity_array():
            all_i_obs.append(array)
        if labels == "BATCH":
            assert batches is None
            batches = array
    if i_obs is None:
        if len(all_i_obs) == 0:
            raise Sorry("No intensities found")
        elif len(all_i_obs) > 1:
            if labels is not None:
                lab_tab = label_table(all_i_obs)
                i_obs = lab_tab.select_array(
                    label=labels[0], command_line_switch="labels"
                )
            if i_obs is None:
                raise Sorry(
                    "Multiple intensity arrays - please specify one:\n%s"
                    % "\n".join(
                        ["  labels=%s" % a.info().label_string() for a in all_i_obs]
                    )
                )
        else:
            i_obs = all_i_obs[0]
    # need original miller indices otherwise we don't get correct anomalous
    # merging statistics
    if "M_ISYM" in mtz_object.column_labels():
        indices = mtz_object.extract_original_index_miller_indices()
        i_obs = i_obs.customized_copy(indices=indices, info=i_obs.info())
    return i_obs, batches


phil_str = """
  rmerge = None
    .type = float(value_min=0)
    .help = "Maximum value of Rmerge in the outer resolution shell"
    .short_caption = "Outer shell Rmerge"
    .expert_level = 1
  completeness = None
    .type = float(value_min=0)
    .help = "Minimum completeness in the outer resolution shell"
    .short_caption = "Outer shell completeness"
    .expert_level = 1
  cc_ref = 0.1
    .type = float(value_min=0)
    .help = "Minimum value of CC vs reference data set in the outer resolution shell"
    .short_caption = "Outer shell CCref"
    .expert_level = 1
  cc_half = 0.3
    .type = float(value_min=0)
    .help = "Minimum value of CC½ in the outer resolution shell"
    .short_caption = "Outer shell CC½"
    .expert_level = 1
  cc_half_method = *half_dataset sigma_tau
    .type = choice
    .short_caption = "CC½ method"
  cc_half_significance_level = 0.1
    .type = float(value_min=0, value_max=1)
    .expert_level = 1
    .short_caption = "CC½ significance level"
  cc_half_fit = polynomial *tanh
    .type = choice
    .expert_level = 1
    .short_caption = "CC½ fit"
  isigma = None
    .type = float(value_min=0)
    .help = "Minimum value of the unmerged <I/sigI> in the outer resolution shell"
    .short_caption = "Outer shell unmerged <I/sigI>"
    .expert_level = 1
  misigma = None
    .type = float(value_min=0)
    .help = "Minimum value of the merged <I/sigI> in the outer resolution shell"
    .short_caption = "Outer shell merged <I/sigI>"
    .expert_level = 1
  i_mean_over_sigma_mean = None
    .type = float(value_min=0)
    .help = "Minimum value of the unmerged <I>/<sigI> in the outer resolution shell"
    .short_caption = "Outer shell unmerged <I>/<sigI>"
    .expert_level = 2
  nbins = 100
    .type = int
    .help = "Maximum number of resolution bins to use for estimation of resolution limit."
    .short_caption = "Number of resolution bins."
    .expert_level = 1
  reflections_per_bin = 10
    .type = int
    .help = "Minimum number of reflections per bin."
    .short_caption = "Minimum number of reflections per bin"
  binning_method = *counting_sorted volume
    .type = choice
    .help = "Use equal-volume bins or bins with approximately equal numbers of reflections per bin."
    .short_caption = "Equal-volume or equal #ref binning."
    .expert_level = 1
  anomalous = False
    .type = bool
    .help = "Keep anomalous pairs separate in merging statistics"
    .short_caption = "Anomalous"
    .expert_level = 1
  labels = None
    .type = strings
    .short_caption = "Labels"
  space_group = None
    .type = space_group
    .expert_level = 1
    .short_caption = "Space group"
  reference = None
    .type = path
    .short_caption = "Reference"
  emax = 4
    .type = float(value_min = 0)
    .help = "Reject reflections with normalised intensities E^2 > emax^2"
    .short_caption = "Maximum normalised intensity"
"""


phil_defaults = iotbx.phil.parse(
    """
resolution {
%s
  batch_range = None
    .type = ints(size=2, value_min=0)
}
"""
    % phil_str
)


def plot_result(metric, result):
    if metric == metrics.CC_HALF:
        return plots.cc_half_plot(
            result.d_star_sq,
            result.y_obs,
            cc_half_critical_values=result.critical_values,
            cc_half_fit=result.y_fit,
            d_min=result.d_min,
        )
    else:
        d = {
            metrics.MISIGMA: "Merged <I/σ(I)>",
            metrics.ISIGMA: "Unmerged <I/σ(I)>",
            metrics.I_MEAN_OVER_SIGMA_MEAN: "&lt;I&gt;/<σ(I)>",
            metrics.RMERGE: "R<sub>merge</sub> ",
            metrics.COMPLETENESS: "Completeness",
        }
        d_star_sq_tickvals, d_star_sq_ticktext = plots.d_star_sq_to_d_ticks(
            result.d_star_sq, 5
        )
        return {
            "data": [
                {
                    "x": list(result.d_star_sq),  # d_star_sq
                    "y": list(result.y_obs),
                    "type": "scatter",
                    "name": "y_obs",
                },
                (
                    {
                        "x": list(result.d_star_sq),
                        "y": list(result.y_fit),
                        "type": "scatter",
                        "name": "y_fit",
                        "line": {"color": "rgb(47, 79, 79)"},
                    }
                    if result.y_fit
                    else {}
                ),
                (
                    {
                        "x": [uctbx.d_as_d_star_sq(result.d_min)] * 2,
                        "y": [
                            0,
                            max(
                                1,
                                flex.max(result.y_obs),
                                flex.max(result.y_fit) if result.y_fit else 0,
                            ),
                        ],
                        "type": "scatter",
                        "name": f"d_min = {result.d_min:.2f} Å",
                        "mode": "lines",
                        "line": {"color": "rgb(169, 169, 169)", "dash": "dot"},
                    }
                    if result.d_min
                    else {}
                ),
            ],
            "layout": {
                "title": f"{d.get(metric)} vs. resolution",
                "xaxis": {
                    "title": "Resolution (Å)",
                    "tickvals": d_star_sq_tickvals,
                    "ticktext": d_star_sq_ticktext,
                },
                "yaxis": {"title": d.get(metric), "rangemode": "tozero"},
            },
        }


class ResolutionResult(typing.NamedTuple):
    d_star_sq: flex.double
    y_obs: flex.double
    y_fit: flex.double
    d_min: float
    critical_values: flex.double = None


class Resolutionizer:
    """A class to calculate things from merging reflections."""

    def __init__(self, i_obs, params, batches=None, reference=None):

        self._params = params
        self._reference = reference

        if self._reference is not None:
            self._reference = self._reference.merge_equivalents(
                use_internal_variance=False
            ).array()

        i_obs = i_obs.customized_copy(
            anomalous_flag=params.anomalous, info=i_obs.info()
        )

        if self._params.batch_range is not None and batches is not None:
            batch_min, batch_max = self._params.batch_range
            assert batches is not None
            sel = (batches.data() >= batch_min) & (batches.data() <= batch_max)
            i_obs = i_obs.select(sel).set_info(i_obs.info())

        if self._params.space_group is not None:
            i_obs = i_obs.customized_copy(
                space_group_info=self._params.space_group, info=i_obs.info()
            )

        if self._params.emax:
            normalised = quasi_normalisation(i_obs)
            e2_cutoff = self._params.emax**2
            sel = normalised.data() < e2_cutoff
            logger.info(
                f"Removing {sel.count(False)} Wilson outliers with E^2 >= {e2_cutoff}"
            )
            i_obs = i_obs.select(sel)

        self._intensities = i_obs

        self._merging_statistics = iotbx.merging_statistics.dataset_statistics(
            i_obs=i_obs,
            n_bins=self._params.nbins,
            reflections_per_bin=self._params.reflections_per_bin,
            cc_one_half_significance_level=self._params.cc_half_significance_level,
            cc_one_half_method=self._params.cc_half_method,
            binning_method=self._params.binning_method,
            anomalous=params.anomalous,
            use_internal_variance=False,
            eliminate_sys_absent=False,
            assert_is_not_unique_set_under_symmetry=False,
        )

    @classmethod
    def from_unmerged_mtz(cls, scaled_unmerged, params):
        """Construct the resolutionizer from an mtz file."""

        i_obs, batches = miller_array_from_mtz(
            scaled_unmerged, anomalous=params.anomalous, labels=params.labels
        )
        if params.reference is not None:
            reference, _ = miller_array_from_mtz(
                params.reference, anomalous=params.anomalous, labels=params.labels
            )
        else:
            reference = None

        return cls(i_obs, params, batches=batches, reference=reference)

    @classmethod
    def from_reflections_and_experiments(cls, reflection_tables, experiments, params):
        """Construct the resolutionizer from native dials datatypes."""
        # add some assertions about data

        # do batch assignment (same functions as in dials.export)
        offsets = calculate_batch_offsets(experiments)
        reflection_tables = assign_batches_to_reflections(reflection_tables, offsets)
        batches = flex.int()
        intensities = flex.double()
        indices = flex.miller_index()
        variances = flex.double()
        for table in reflection_tables:
            if "intensity.scale.value" in table:
                table = filter_reflection_table(
                    table, ["scale"], partiality_threshold=0.4
                )
                intensities.extend(table["intensity.scale.value"])
                variances.extend(table["intensity.scale.variance"])
            else:
                table = filter_reflection_table(
                    table, ["profile"], partiality_threshold=0.4
                )
                intensities.extend(table["intensity.prf.value"])
                variances.extend(table["intensity.prf.variance"])
            indices.extend(table["miller_index"])
            batches.extend(table["batch"])

        crystal_symmetry = miller.crystal.symmetry(
            unit_cell=determine_best_unit_cell(experiments),
            space_group=experiments[0].crystal.get_space_group(),
            assert_is_compatible_unit_cell=False,
        )
        miller_set = miller.set(crystal_symmetry, indices, anomalous_flag=False)
        i_obs = miller.array(miller_set, data=intensities, sigmas=flex.sqrt(variances))
        i_obs.set_observation_type_xray_intensity()
        i_obs.set_info(miller.array_info(source="DIALS", source_type="refl"))

        ms = i_obs.customized_copy()
        batch_array = miller.array(ms, data=batches)

        if params.reference is not None:
            reference, _ = miller_array_from_mtz(
                params.reference, anomalous=params.anomalous, labels=params.labels
            )
        else:
            reference = None

        return cls(i_obs, params, batches=batch_array, reference=reference)

    def resolution(self, metric, limit=None):
        if metric == metrics.CC_HALF:
            return resolution_cc_half(
                self._merging_statistics,
                limit,
                cc_half_method=self._params.cc_half_method,
                model=tanh_fit
                if self._params.cc_half_fit == "tanh"
                else polynomial_fit,
            )
        elif metric == metrics.CC_REF:
            return self._resolution_cc_ref(limit=self._params.cc_ref)
        else:
            model = {
                metrics.RMERGE: log_inv_fit,
                metrics.COMPLETENESS: polynomial_fit,
                metrics.ISIGMA: log_fit,
                metrics.MISIGMA: log_fit,
                metrics.I_MEAN_OVER_SIGMA_MEAN: log_fit,
            }[metric]
            return resolution_fit_from_merging_stats(
                self._merging_statistics, metric.value, model, limit
            )

    def resolution_auto(self):
        """Compute resolution limits based on the current self._params set."""

        metric_to_output = {
            metrics.ISIGMA: "I/sig",
            metrics.MISIGMA: "Mn(I/sig)",
            metrics.I_MEAN_OVER_SIGMA_MEAN: "Mn(I)/Mn(sig)",
        }

        plot_d = {}

        for metric in metrics:
            name = metric.name.lower()
            limit = getattr(self._params, name)
            if metric == metrics.CC_REF and not self._reference:
                limit = None
            if limit:
                try:
                    result = self.resolution(metric, limit=limit)
                except RuntimeError as e:
                    logger.info(f"Resolution fit against {name} failed: {e}")
                    continue
                pretty_name = metric_to_output.get(metric, name)
                if result.d_min:
                    logger.info(
                        f"Resolution {pretty_name}:{result.d_min:{18 - len(pretty_name)}.2f}"
                    )
                plot_d[name] = plot_result(metric, result)
        return plot_d

    def _resolution_cc_ref(self, limit=None):
        """Compute a resolution limit where cc_ref < 0.5 (limit if
        set) or the full extent of the data."""

        if limit is None:
            limit = self._params.cc_ref

        intensities = self._intensities.merge_equivalents(
            use_internal_variance=False
        ).array()
        cc_s = flex.double()
        for b in self._merging_statistics.bins:
            cc = intensities.resolution_filter(
                d_min=b.d_min, d_max=b.d_max
            ).correlation(
                self._reference.resolution_filter(d_min=b.d_min, d_max=b.d_max),
                assert_is_similar_symmetry=False,
            )
            cc_s.append(cc.coefficient())
        cc_s = cc_s.reversed()

        fit = tanh_fit if self._params.cc_half_fit == "tanh" else polynomial_fit
        d_star_sq = flex.double(
            1 / b.d_min**2 for b in self._merging_statistics.bins
        ).reversed()

        return resolution_fit(d_star_sq, cc_s, fit, limit)
