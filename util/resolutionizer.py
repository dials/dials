from __future__ import absolute_import, division, print_function

import logging
import math
from dataclasses import dataclass

import iotbx.mtz
import iotbx.phil
import iotbx.merging_statistics
from cctbx.array_family import flex
from cctbx import miller, uctbx
from iotbx.reflection_file_utils import label_table
from scitbx.math import curve_fitting
from scitbx.math import five_number_summary

from dials.algorithms.scaling.scaling_library import determine_best_unit_cell
from dials.util import Sorry
from dials.util.batch_handling import (
    calculate_batch_offsets,
    assign_batches_to_reflections,
)
from dials.util.filter_reflections import filter_reflection_table
from dials.util import tabulate

logger = logging.getLogger(__name__)


def polynomial_fit(x, y, degree=5):
    """Fit the values y(x) then return this fit. x, y should
    be iterables containing floats of the same size. The order is the order
    of polynomial to use for this fit. This will be useful for e.g. I/sigma."""

    fit = curve_fitting.univariate_polynomial_fit(x, flex.log(y), degree=degree)
    f = curve_fitting.univariate_polynomial(*fit.params)
    return flex.exp(f(x))


def tanh_fit(x, y, iqr_multiplier=None):
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
    """Fit the values log(y(x)) then return exp() to this fit. x, y should
    be iterables containing floats of the same size. The order is the order
    of polynomial to use for this fit. This will be useful for e.g. I/sigma."""

    fit = curve_fitting.univariate_polynomial_fit(x, flex.log(y), degree=degree)
    f = curve_fitting.univariate_polynomial(*fit.params)
    return flex.exp(f(x))


def log_inv_fit(x, y, degree=5):
    """Fit the values log(1 / y(x)) then return the inverse of this fit.
    x, y should be iterables, the order of the polynomial for the transformed
    fit needs to be specified. This will be useful for e.g. Rmerge."""

    fit = curve_fitting.univariate_polynomial_fit(x, flex.log(1 / y), degree=degree)
    f = curve_fitting.univariate_polynomial(*fit.params)
    return 1 / flex.exp(f(x))


def resolution_fit_from_merging_stats(merging_stats, metric, model, limit, sel=None):
    y_obs = flex.double(getattr(b, metric) for b in merging_stats.bins).reversed()
    d_star_sq = flex.double(
        uctbx.d_as_d_star_sq(b.d_min) for b in merging_stats.bins
    ).reversed()
    return resolution_fit(d_star_sq, y_obs, model, limit, sel=sel)


def resolution_fit(d_star_sq, y_obs, model, limit, sel=None):
    if not sel:
        sel = flex.bool(len(d_star_sq), True)
    sel &= y_obs > 0
    y_obs = y_obs.select(sel)
    d_star_sq = d_star_sq.select(sel)

    if flex.min(y_obs) > limit:
        d_min = 1.0 / math.sqrt(flex.max(d_star_sq))
        y_fit = None

    else:
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

        # rlimit = limit * max(y_obs)
        try:
            d_min = 1.0 / math.sqrt(interpolate_value(d_star_sq, y_fit, limit))
        except Exception as e:
            logger.info(e)
            d_min = uctbx.d_star_sq_as_d(flex.max(d_star_sq))

    return ResolutionResult(d_star_sq, y_obs, y_fit, d_min)


def get_cc_half_significance(merging_stats, cc_half_method):
    if cc_half_method == "sigma_tau":
        return flex.bool(
            b.cc_one_half_sigma_tau_significance for b in merging_stats.bins
        ).reversed()
    else:
        return flex.bool(
            b.cc_one_half_significance for b in merging_stats.bins
        ).reversed()


def get_cc_half_critical_values(merging_stats, cc_half_method):
    if cc_half_method == "sigma_tau":
        return flex.double(
            b.cc_one_half_sigma_tau_critical_value for b in merging_stats.bins
        ).reversed()
    else:
        return flex.double(
            b.cc_one_half_critical_value for b in merging_stats.bins
        ).reversed()


def resolution_cc_half(
    merging_stats, limit, cc_half_method="half_dataset", model=tanh_fit
):
    """Compute a resolution limit where cc_half < 0.5 (limit if
    set) or the full extent of the data."""

    sel = get_cc_half_significance(merging_stats, cc_half_method)
    metric = "cc_one_half_sigma_tau" if cc_half_method == "sigma_tau" else "cc_one_half"
    result = resolution_fit_from_merging_stats(
        merging_stats, metric, model, limit, sel=sel
    )
    result.critical_values = get_cc_half_critical_values(
        merging_stats, cc_half_method
    ).select(sel)
    return result


def interpolate_value(x, y, t):
    """Find the value of x: y(x) = t."""

    if t > max(y) or t < min(y):
        raise RuntimeError("t outside of [%f, %f]" % (min(y), max(y)))

    for j in range(1, len(x)):
        x0 = x[j - 1]
        y0 = y[j - 1]

        x1 = x[j]
        y1 = y[j]

        if (y0 - t) * (y1 - t) < 0:
            return x0 + (t - y0) * (x1 - x0) / (y1 - y0)


def miller_array_from_mtz(unmerged_mtz, params):
    mtz_object = iotbx.mtz.object(file_name=unmerged_mtz)
    miller_arrays = mtz_object.as_miller_arrays(
        merge_equivalents=False, anomalous=params.anomalous
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
            if params.labels is not None:
                lab_tab = label_table(all_i_obs)
                i_obs = lab_tab.select_array(
                    label=params.labels[0], command_line_switch="labels"
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
    .help = "Minimum value of CC vs reference dataset in the outer resolution shell"
    .short_caption = "Outer shell CCref"
    .expert_level = 1
  cc_half = 0.3
    .type = float(value_min=0)
    .help = "Minimum value of CC1/2 in the outer resolution shell"
    .short_caption = "Outer shell CC1/2"
    .expert_level = 1
  cc_half_method = *half_dataset sigma_tau
    .type = choice
  cc_half_significance_level = 0.1
    .type = float(value_min=0, value_max=1)
    .expert_level = 1
  cc_half_fit = polynomial *tanh
    .type = choice
    .expert_level = 1
  isigma = 0.25
    .type = float(value_min=0)
    .help = "Minimum value of the unmerged <I/sigI> in the outer resolution shell"
    .short_caption = "Outer shell unmerged <I/sigI>"
    .expert_level = 1
  misigma = 1.0
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
  binning_method = *counting_sorted volume
    .type = choice
    .help = "Use equal-volume bins or bins with approximately equal numbers of reflections per bin."
    .short_caption = "Equal-volume or equal #ref binning."
    .expert_level = 1
  anomalous = False
    .type = bool
    .short_caption = "Keep anomalous pairs separate in merging statistics"
    .expert_level = 1
  labels = None
    .type = strings
  space_group = None
    .type = space_group
    .expert_level = 1
  reference = None
    .type = path
"""


phil_defaults = iotbx.phil.parse(
    """
resolutionizer {
%s
  batch_range = None
    .type = ints(size=2, value_min=0)
  plot = False
    .type = bool
    .expert_level = 2
}
"""
    % phil_str
)


def plot_resolution_result(result, ylabel, ax=None):
    from matplotlib import pyplot as plt

    plt.style.use("ggplot")
    ylabel = ylabel
    if ax is None:
        fig, ax = plt.subplots()
    else:
        ax = ax.figure

    if result.y_fit:
        ax.plot(result.d_star_sq, result.y_fit, label="Fit")
    ax.plot(result.d_star_sq, result.y_obs, label="Raw")
    if result.critical_values:
        ax.plot(result.d_star_sq, result.critical_values, label="Critical values")
    if ylabel.startswith("CC"):
        ylim = ax.get_ylim()
        ax.set_ylim(0, max(ylim[1], 1.05))

    d_star_sq = uctbx.d_as_d_star_sq(result.d_min)
    ax.plot([d_star_sq, d_star_sq], ax.get_ylim(), linestyle="--")

    xticks = ax.get_xticks()
    xticks_d = ["%.2f" % uctbx.d_star_sq_as_d(ds2) if ds2 > 0 else 0 for ds2 in xticks]
    ax.set_xticklabels(xticks_d)
    ax.set_xlabel("Resolution (A)")
    ax.set_ylabel(ylabel)
    ax.legend(loc="best")
    return ax


@dataclass
class ResolutionResult:
    d_star_sq: flex.double
    y_obs: flex.double
    y_fit: flex.double
    d_min: float
    critical_values: flex.double = None


class Resolutionizer(object):
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

        i_obs, batches = miller_array_from_mtz(scaled_unmerged, params)
        if params.reference is not None:
            reference, _ = miller_array_from_mtz(params.reference, params)
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
            reference, _ = miller_array_from_mtz(params.reference, params)
        else:
            reference = None

        return cls(i_obs, params, batches=batch_array, reference=reference)

    def resolution_auto(self):
        """Compute resolution limits based on the current self._params set."""

        if self._params.plot:
            import matplotlib

            matplotlib.use("Agg")

            from matplotlib import pyplot as plt

        if self._params.rmerge:
            result = self.resolution_rmerge()
            if self._params.plot:
                plot_resolution_result(result, "Rmerge")
                plt.savefig("rmerge.png")
            logger.info("Resolution rmerge:       %.2f", result.d_min)

        if self._params.completeness:
            result = self.resolution_completeness()
            if self._params.plot:
                plot_resolution_result(result, "Completeness")
                plt.savefig("completeness.png")
            logger.info("Resolution completeness: %.2f", result.d_min)

        if self._params.cc_half:
            result = self.resolution_cc_half()
            if self._params.plot:
                plot_resolution_result(result, "CC½")
                plt.savefig("cc_half.png")
            logger.info("Resolution cc_half:      %.2f", result.d_min)

        if self._params.cc_ref and self._reference is not None:
            result = self.resolution_cc_ref()
            if self._params.plot:
                plot_resolution_result(result, "CCref")
                plt.savefig("cc_ref.png")
            logger.info("Resolution cc_ref:       %.2f", result.d_min)

        if self._params.isigma:
            result = self.resolution_unmerged_isigma()
            if self._params.plot:
                plot_resolution_result(result, "Unmerged I/sigma")
                plt.savefig("isigma.png")
            logger.info("Resolution I/sig:        %.2f", result.d_min)

        if self._params.misigma:
            result = self.resolution_merged_isigma()
            if self._params.plot:
                plot_resolution_result(result, "Merged I/sigma")
                plt.savefig("misigma.png")
            logger.info("Resolution Mn(I/sig):    %.2f", result.d_min)

        if self._params.i_mean_over_sigma_mean:
            result = self.resolution_i_mean_over_sigma_mean()
            if self._params.plot:
                plot_resolution_result(result, "Unmerged <I>/<sigma>")
                plt.savefig("i_mean_over_sigma_mean.png")
            logger.info(
                "Resolution Mn(I)/Mn(sig):    %.2f", result.d_min,
            )

    def resolution_rmerge(self, limit=None):
        """Compute a resolution limit where either rmerge = 1.0 (limit if
        set) or the full extent of the data. N.B. this fit is only meaningful
        for positive values."""

        if limit is None:
            limit = self._params.rmerge

        return resolution_fit_from_merging_stats(
            self._merging_statistics, "r_merge", log_inv_fit, limit
        )

    def resolution_i_mean_over_sigma_mean(self, limit=None):
        """Compute a resolution limit where either <I>/<sigma> = 1.0 (limit if
        set) or the full extent of the data."""

        if limit is None:
            limit = self._params.i_mean_over_sigma_mean

        return resolution_fit_from_merging_stats(
            self._merging_statistics, "i_mean_over_sigi_mean", log_fit, limit
        )

    def resolution_unmerged_isigma(self, limit=None):
        """Compute a resolution limit where either I/sigma = 1.0 (limit if
        set) or the full extent of the data."""

        if limit is None:
            limit = self._params.isigma

        return resolution_fit_from_merging_stats(
            self._merging_statistics, "unmerged_i_over_sigma_mean", log_fit, limit
        )

    def resolution_merged_isigma(self, limit=None):
        """Compute a resolution limit where either Mn(I/sigma) = 1.0 (limit if
        set) or the full extent of the data."""

        if limit is None:
            limit = self._params.misigma

        return resolution_fit_from_merging_stats(
            self._merging_statistics, "i_over_sigma_mean", log_fit, limit
        )

    def resolution_completeness(self, limit=None):
        """Compute a resolution limit where completeness < 0.5 (limit if
        set) or the full extent of the data. N.B. this completeness is
        with respect to the *maximum* completeness in a shell, to reflect
        triclinic cases."""

        if limit is None:
            limit = self._params.completeness

        return resolution_fit_from_merging_stats(
            self._merging_statistics, "completeness", polynomial_fit, limit
        )

    def resolution_cc_half(self, limit=None):
        """Compute a resolution limit where cc_half < 0.5 (limit if
        set) or the full extent of the data."""

        if limit is None:
            limit = self._params.cc_half
        return resolution_cc_half(
            self._merging_statistics,
            limit,
            cc_half_method=self._params.cc_half_method,
            model=tanh_fit if self._params.cc_half_fit == "tanh" else polynomial_fit,
        )

    def resolution_cc_ref(self, limit=None):
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
            1 / b.d_min ** 2 for b in self._merging_statistics.bins
        ).reversed()

        return resolution_fit(d_star_sq, cc_s, fit, limit)
