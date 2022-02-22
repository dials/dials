"""
Error model classes for scaling.
"""

from __future__ import annotations

import logging
from math import ceil, exp, log

import numpy as np
from scipy.stats import norm

from dxtbx import flumpy
from iotbx import phil
from scitbx import sparse

from dials.array_family import flex
from dials.util import tabulate

logger = logging.getLogger("dials")

phil_scope = phil.parse(
    """
    error_model = *basic None
        .type = choice
        .help = "The error model to use."
        .expert_level = 1
    basic {
        a = None
          .type = float
          .help = "Used this fixed value for the error model 'a' parameter"
          .expert_level = 2
        b = None
          .type = float
          .help = "Used this fixed value for the error model 'b' parameter"
          .expert_level = 2
        minimisation = *individual regression None
          .type = choice
          .help = "The algorithm to use for basic error model minimisation."
                  "For individual, the a and b parameters are optimised"
                  "sequentially. For regression, a linear fit is made to"
                  "determine both parameters concurrently. If minimisation=None,"
                  "the model parameters are fixed to their initial or given values."
          .expert_level = 3
        min_Ih = 25.0
            .type = float
            .help = "Reflections with expected intensity above this value are to."
                    "be used in error model minimisation."
            .expert_level = 2
        n_bins = 10
            .type = int
            .help = "The number of intensity bins to use for the error model optimisation."
            .expert_level = 2
    }
    reset_error_model = False
        .type = bool
        .help = "If True, the error model is reset to the default at the start"
                "of scaling, as opposed to loading the current error model."
    grouping = individual grouped *combined
        .type = choice
        .help = "This options selects whether one error model is determined"
                "for all sweeps (combined), whether one error model is"
                "determined per-sweep (individual), or whether a custom"
                "grouping should be used. If grouping=grouped, each group"
                "should be specified with the error_model_group=parameter."
        .expert_level = 2
    error_model_group = None
        .type = ints
        .multiple = True
        .help = "Specify a subset of sweeps which should share an error model."
                "If no groups are specified here, this is interpreted to mean"
                "that all sweeps should share a common error model."

    """
)


def calc_sigmaprime(x, Ih_table) -> np.array:
    """Calculate the error from the model."""
    sigmaprime = (
        x[0] * np.sqrt(Ih_table.variances + np.square(x[1] * Ih_table.intensities))
    ) / Ih_table.inverse_scale_factors
    return sigmaprime


def calc_deltahl(Ih_table, n_h, sigmaprime) -> np.array:
    """Calculate the normalised deviations from the model."""
    I_hl = Ih_table.intensities
    g_hl = Ih_table.inverse_scale_factors
    I_h = Ih_table.Ih_values
    prefactor = np.sqrt((n_h - np.full(n_h.size, 1.0)) / n_h)
    delta_hl = prefactor * ((I_hl / g_hl) - I_h) / sigmaprime
    return delta_hl


class ErrorModelRegressionAPM:

    """Parameter manager for error model minimisation using the linear
    regression method.

    Allows refining of just a, just b or both."""

    def __init__(self, model, active_parameters):
        self.model = model
        self.components = {}
        self.active_parameters = active_parameters  # e.g. ["a", "b"]
        self.x = flex.double([])
        n_cumul_params = 0
        for p in active_parameters:
            if p in model.components:
                self.x.extend(model.components[p].parameters)
                self.components.update(
                    {p: {"start_idx": n_cumul_params, "end_idx": len(self.x)}}
                )
                n_cumul_params += 1

    def select_parameters(self, component):
        """Select the subset of self.x corresponding to the component (a string)."""
        start_idx = self.components[component]["start_idx"]
        end_idx = self.components[component]["end_idx"]
        return self.x[start_idx:end_idx]

    def set_param_vals(self, x):
        """Set method for refinement engine access."""
        self.x = x
        for name, component in self.model.components.items():
            if name in self.active_parameters:
                component.parameters = self.select_parameters(name)

    def get_param_vals(self):
        """Get method for refinement engine access."""
        return self.x

    def resolve_model_parameters(self):
        """Update parameters in the model."""
        if self.active_parameters == ["a"]:
            if not self.x[0] > 0:
                raise ValueError(
                    "Error model refinement resulted in a negative value for a^2"
                )
            a = self.x[0] ** 0.5
            b = self.model.parameters[1]
        elif self.active_parameters == ["b"]:
            if not self.x[0] > 0:
                raise ValueError(
                    "Error model refinement resulted in a negative value for b^2"
                )
            a = self.model.parameters[0]
            b = self.x[0] ** 0.5
        else:
            if not self.x[0] > 0:
                raise ValueError(
                    "Error model refinement resulted in a negative value for a^2"
                )
            if not self.x[1] > 0:
                raise ValueError(
                    "Error model refinement resulted in a negative value for b^2"
                )
            a = self.x[0] ** 0.5
            b = self.x[1] ** 0.5 / a
        self.model.update([a, b])


class ErrorModelA_APM:

    """Parameter manager for minimising A component with individual minimizer"""

    def __init__(self, model):
        self.model = model
        self.x = [0.0, model.components["a"].parameters[0]]

    def set_param_vals(self, x):
        """method for refinement engine access"""
        self.x = x

    def get_param_vals(self):
        """method for refinement engine access"""
        return self.x

    def resolve_model_parameters(self):
        """Update parameters in the model."""
        self.model.components["a"].parameters *= self.x[1]


class ErrorModelB_APM:

    """Parameter manager for minimising Bcomponent with individual minimizer"""

    def __init__(self, model):
        self.model = model
        self.x = [model.components["b"].parameters[0]]

    def set_param_vals(self, x):
        """method for refinement engine access"""
        self.x = x
        self.model.components["b"].parameters = flex.double([self.x[0]])
        self.model.binner.update(self.model.parameters)

    def get_param_vals(self):
        """method for refinement engine access"""
        return self.x

    def resolve_model_parameters(self):
        """Update parameters in the model."""
        self.model.components["b"].parameters = flex.double([self.x[0]])


class ErrorModelBinner:

    """A binner for the error model data.

    Data are binned based on Ih, and methods are available for
    calculating bin variances, summation within bins etc."""

    def __init__(self, Ih_table, min_reflections_required=250, n_bins=10):
        self.binning_info = {
            "initial_variances": [],
            "bin_boundaries": [],
            "mean_intensities": np.array([], dtype=float).reshape((0,)),
            "bin_variances": [],
            "refl_per_bin": [],
            "n_reflections": None,
        }
        self.n_bins = n_bins
        self.Ih_table = Ih_table
        self.min_reflections_required = min_reflections_required
        self.n_h = self.Ih_table.calc_nh()
        self.sigmaprime = calc_sigmaprime([1.0, 0.0], self.Ih_table)
        self.summation_matrix = self._create_summation_matrix()
        self.weights = np.array(self.binning_info["mean_intensities"])
        self.delta_hl = calc_deltahl(self.Ih_table, self.n_h, self.sigmaprime)
        self.bin_variances = self.calculate_bin_variances()
        self.binning_info["initial_variances"] = self.binning_info["bin_variances"]

    def update(self, parameters):
        """Update the variances for updated model parameters."""
        self.sigmaprime = calc_sigmaprime(parameters, self.Ih_table)
        self.delta_hl = calc_deltahl(self.Ih_table, self.n_h, self.sigmaprime)
        self.bin_variances = self.calculate_bin_variances()

    def _create_summation_matrix(self):
        """Create a summation matrix to allow sums into intensity bins.

        This routine attempts to bin into bins equally spaced in log(intensity),
        to give a representative sample across all intensities. To avoid
        undersampling, it is required that there are at least 100 reflections
        per intensity bin unless there are very few reflections."""
        n = self.Ih_table.size
        self.binning_info["n_reflections"] = n
        summation_matrix = sparse.matrix(n, self.n_bins)
        Ih = self.Ih_table.Ih_values * self.Ih_table.inverse_scale_factors
        size_order = flex.sort_permutation(flumpy.from_numpy(Ih), reverse=True)
        Imax = Ih.max()
        min_Ih = Ih.min()
        Imin = max(1.0, min_Ih)  # avoid log issues
        spacing = (log(Imax) - log(Imin)) / float(self.n_bins)
        boundaries = [Imax] + [
            exp(log(Imax) - (i * spacing)) for i in range(1, self.n_bins + 1)
        ]
        boundaries[-1] = min_Ih - 0.01
        self.binning_info["bin_boundaries"] = np.array(boundaries)
        self.binning_info["refl_per_bin"] = np.array([], dtype=float).reshape((0,))

        n_cumul = 0
        if Ih.size > 100 * self.min_reflections_required:
            self.min_reflections_required = int(Ih.size / 100.0)
        min_per_bin = min(self.min_reflections_required, int(n / (3.0 * self.n_bins)))
        for i in range(len(self.binning_info["bin_boundaries"]) - 1):
            maximum = self.binning_info["bin_boundaries"][i]
            minimum = self.binning_info["bin_boundaries"][i + 1]
            sel1 = Ih <= maximum
            sel2 = Ih > minimum
            sel = sel1 & sel2
            # isel = sel.iselection()
            n_in_bin = np.count_nonzero(sel)  # isel.size()
            if n_in_bin < min_per_bin:  # need more in this bin
                m = n_cumul + min_per_bin
                if m < n:  # still some refl left to use
                    idx = size_order[m]
                    intensity = Ih[idx]
                    self.binning_info["bin_boundaries"][i + 1] = intensity
                    minimum = self.binning_info["bin_boundaries"][i + 1]
                    sel = sel1 & (Ih > minimum)
                    # isel = sel.iselection()
                    n_in_bin = np.count_nonzero(sel)  # isel.size()
            self.binning_info["refl_per_bin"] = np.append(
                self.binning_info["refl_per_bin"], [n_in_bin]
            )
            for j in np.nonzero(sel)[0]:  # isel:
                summation_matrix[int(j), i] = 1
            n_cumul += n_in_bin
        cols_to_del = []
        for i, col in enumerate(summation_matrix.cols()):
            if col.non_zeroes < min_per_bin - 5:
                cols_to_del.append(i)
        n_new_cols = summation_matrix.n_cols - len(cols_to_del)
        if n_new_cols == self.n_bins:
            for i in range(len(self.binning_info["bin_boundaries"]) - 1):
                maximum = self.binning_info["bin_boundaries"][i]
                minimum = self.binning_info["bin_boundaries"][i + 1]
                sel1 = Ih <= maximum
                sel2 = Ih > minimum
                sel = sel1 & sel2
                self.binning_info["mean_intensities"] = np.append(
                    self.binning_info["mean_intensities"], [np.mean(Ih[sel])]
                )
            return summation_matrix
        new_sum_matrix = sparse.matrix(summation_matrix.n_rows, n_new_cols)
        next_col = 0
        refl_per_bin = np.array([], dtype=float).reshape((0,))
        new_bounds = np.array([], dtype=float).reshape((0,))
        for i, col in enumerate(summation_matrix.cols()):
            if i not in cols_to_del:
                new_sum_matrix[:, next_col] = col
                next_col += 1
                new_bounds = np.append(
                    new_bounds, [self.binning_info["bin_boundaries"][i]]
                )
                refl_per_bin = np.append(
                    refl_per_bin, [self.binning_info["refl_per_bin"][i]]
                )
        self.binning_info["refl_per_bin"] = refl_per_bin
        new_bounds = np.append(new_bounds, [self.binning_info["bin_boundaries"][-1]])
        self.binning_info["bin_boundaries"] = new_bounds
        for i in range(len(new_bounds) - 1):
            maximum = new_bounds[i]
            minimum = new_bounds[i + 1]
            sel1 = Ih <= maximum
            sel2 = Ih > minimum
            sel = sel1 & sel2
            self.binning_info["mean_intensities"] = np.append(
                self.binning_info["mean_intensities"], [np.mean(Ih[sel])]
            )
        return new_sum_matrix

    def calculate_bin_variances(self) -> np.array:
        """Calculate the variance of each bin."""
        sum_deltasq = flumpy.to_numpy(
            flumpy.from_numpy(np.square(self.delta_hl)) * self.summation_matrix
        )
        sum_delta_sq = np.square(
            flumpy.to_numpy(flumpy.from_numpy(self.delta_hl) * self.summation_matrix)
        )
        bin_vars = (sum_deltasq / self.binning_info["refl_per_bin"]) - (
            sum_delta_sq / np.square(self.binning_info["refl_per_bin"])
        )
        self.binning_info["bin_variances"] = bin_vars
        return bin_vars


class BComponent:

    """The basic error model B parameter component"""

    def __init__(self, initial_value=0.02):
        self.parameters = flex.double([initial_value])
        self._n_params = 1


class AComponent:

    """The basic error model A parameter component"""

    def __init__(self, initial_value=1.00):
        self.parameters = flex.double([initial_value])
        self._n_params = 1


class BasicErrorModel:

    """Definition of a basic two-parameter error model."""

    min_reflections_required = 250

    id_ = "basic"

    def __init__(self, a=None, b=None, basic_params=None):

        """
        A basic two-parameter error model s'^2 = a^2(s^2 + (bI)^2)

        If a and b are not given as arguments, the params scope is checked to
        see if a user specified fixed value is set. If no fixed values are given
        then the model starts with the default parameters a=1.0 b=0.02
        """

        self.free_components = []
        self.sortedy = None
        self.sortedx = None
        self.binner = None
        if not basic_params:
            basic_params = phil_scope.fetch().extract().basic
        self.params = basic_params
        self.filtered_Ih_table = None
        if not a:
            a = basic_params.a
            if not a:
                a = 1.0
        if not b:
            b = basic_params.b
            if not b:
                b = 0.02
        self.components = {"a": AComponent(a), "b": BComponent(b)}
        self._active_parameters = []
        # if the parameters have been set in the phil scope, then these are to be fixed
        if not basic_params.a:
            self._active_parameters.append("a")
        if not basic_params.b:
            self._active_parameters.append("b")

    def configure_for_refinement(self, Ih_table, min_partiality=0.4):
        """
        Add data to allow error model refinement.

        Raises: ValueError if insufficient reflections left after filtering.
        """
        self.filtered_Ih_table = self.filter_unsuitable_reflections(
            Ih_table, self.params, min_partiality
        )
        # always want binning info so that can calc for output.
        self.binner = ErrorModelBinner(
            self.filtered_Ih_table, self.min_reflections_required, self.params.n_bins
        )
        # need to calculate sorted deltahl for norm dev plotting (and used by
        # individual a-parameter minimiser)
        self.calculate_sorted_deviations(self.parameters)

        self.binner.update(self.parameters)

    @property
    def active_parameters(self):
        return self._active_parameters

    @property
    def parameters(self):
        """A list of the model parameters."""
        return [
            self.components["a"].parameters[0],
            abs(self.components["b"].parameters[0]),
        ]

    @parameters.setter
    def parameters(self, parameters):
        assert len(parameters) == 2
        self.components["a"].parameters = flex.double([parameters[0]])
        self.components["b"].parameters = flex.double([parameters[1]])

    def finalise(self):
        """Perform any actions after minimisation finished."""
        logger.info(self.binned_variances_summary())

    @property
    def n_refl(self):
        """The number of reflections being used in minimisation"""
        return self.filtered_Ih_table.size

    @classmethod
    def filter_unsuitable_reflections(cls, Ih_table, error_params, min_partiality):
        """Filter suitable reflections for minimisation."""
        return filter_unsuitable_reflections(
            Ih_table,
            min_Ih=error_params.min_Ih,
            min_partiality=min_partiality,
            min_reflections_required=cls.min_reflections_required,
        )

    def calculate_sorted_deviations(self, parameters):
        """Sort the x,y data."""
        sigmaprime = calc_sigmaprime(parameters, self.filtered_Ih_table)
        delta_hl = calc_deltahl(
            self.filtered_Ih_table, self.filtered_Ih_table.calc_nh(), sigmaprime
        )
        central_cutoff = 1.5
        n = delta_hl.size
        self.sortedy = np.sort(delta_hl)
        v1 = norm.cdf(-central_cutoff)
        v2 = norm.cdf(central_cutoff)
        idx_cutoff_min = ceil(
            (v1 * n) - 0.5
        )  # first one within the central cutoff range
        idx_cutoff_max = ceil(
            (v2 * n) - 0.5
        )  # first one above the central cutoff range
        central_n = idx_cutoff_max - idx_cutoff_min
        v = np.linspace(
            start=(idx_cutoff_min + 0.5) / n,
            stop=(idx_cutoff_max + 0.5) / n,
            endpoint=False,
            num=central_n,
        )

        self.sortedx = flumpy.from_numpy(norm.ppf(v))
        self.sortedy = flumpy.from_numpy(self.sortedy[idx_cutoff_min:idx_cutoff_max])

    def update(self, parameters):
        """Update the model with new parameters."""
        self.parameters = parameters
        self.binner.update(parameters)
        self.calculate_sorted_deviations(parameters)

    def update_variances(self, variances, intensities):
        """Use the error model parameter to calculate new values for the variances."""
        new_variance = (self.parameters[0] ** 2) * (
            variances + np.square(self.parameters[1] * intensities)
        )
        return new_variance

    def clear_Ih_table(self):
        """Delete the Ih_table, to free memory."""
        if self.binner:
            self.binner.Ih_table = None

    def __str__(self):
        a = abs(self.parameters[0])
        b = abs(self.parameters[1])
        ISa = f"{1.0 / (b * a):.3f}" if (b * a) > 0 else "Unable to estimate"
        return "\n".join(
            (
                "",
                "Error model details:",
                "  Type: basic",
                f"  Parameters: a = {a:.5f}, b = {b:.5f}",
                "  Error model formula: "
                + "\u03C3"
                + "'"
                + "\xb2"
                + " = a"
                + "\xb2"
                + "("
                + "\u03C3\xb2"
                " + (bI)" + "\xb2" + ")",
                "  estimated I/sigma asymptotic limit: %s" % ISa,
                "",
            )
        )

    def binned_variances_summary(self):
        """Generate a summary of the model minimisation for output."""
        header = [
            "Intensity range (<Ih>)",
            "n_refl",
            "Uncorrected variance",
            "Corrected variance",
        ]
        rows = []
        bin_bounds = [f"{i:.2f}" for i in self.binner.binning_info["bin_boundaries"]]
        for i, (initial_var, bin_var, n_refl) in enumerate(
            zip(
                self.binner.binning_info["initial_variances"],
                self.binner.binning_info["bin_variances"],
                self.binner.binning_info["refl_per_bin"],
            )
        ):
            rows.append(
                [
                    bin_bounds[i] + " - " + bin_bounds[i + 1],
                    str(int(n_refl)),
                    str(round(initial_var, 3)),
                    str(round(bin_var, 3)),
                ]
            )
        return "\n".join(
            (
                "Results of error model refinement. Uncorrected and corrected variances",
                "of normalised intensity deviations for given intensity ranges. Variances",
                "are expected to be ~1.0 for reliable errors (sigmas).",
                tabulate(rows, header),
                "",
            )
        )


def filter_unsuitable_reflections(
    Ih_table, min_Ih, min_partiality, min_reflections_required
):
    """
    Choose reflection groups with n_h > 1, as these have deltas of zero by
    definition and will bias the variance calculations. Also, only use groups
    where <Ih> > 25.0, as the assumptions of normally distributed deltas will
    not hold for low <Ih>."""

    if "partiality" in Ih_table.Ih_table:
        sel = Ih_table.Ih_table["partiality"].to_numpy() > min_partiality
        Ih_table = Ih_table.select(sel)

    n = Ih_table.size
    sum_I_over_var = Ih_table.sum_in_groups(Ih_table.intensities / Ih_table.variances)
    n_per_group = Ih_table.sum_in_groups(np.full(n, 1))
    avg_I_over_var = sum_I_over_var / n_per_group
    sel = avg_I_over_var > 0.85
    Ih_table = Ih_table.select_on_groups(sel)
    n_h = Ih_table.calc_nh()
    scaled_Ih = Ih_table.Ih_values * Ih_table.inverse_scale_factors
    # need a scaled min_Ih, where can reasonably expect norm distribution
    # (use min_Ih=25 by default, sigma ~ 5)
    sel2 = scaled_Ih > min_Ih
    # can't calculate a true deviation for groups of 1
    sel3 = n_h > 1.0
    sel4 = Ih_table.intensities > 0.001
    # don't want to include weaker reflections where the background adds
    # significantly to the variances, as these would no longer be normally
    # distributed and skew the fit.
    Ih_table = Ih_table.select(sel2 & sel3 & sel4)
    n = Ih_table.size
    if n < min_reflections_required:
        raise ValueError(
            "Insufficient reflections (%s < %s) to perform error modelling."
            % (n, min_reflections_required)
        )
    n_h = Ih_table.calc_nh()
    # now make sure any left also have n > 1
    sel = n_h > 1.0
    Ih_table = Ih_table.select(sel)

    #  Filter groups with abnormally high internal variances.
    # For a reasonable quality dataset, if b=0.04, a=1.25, then for large Imax,
    # the ratio of the corrected per-reflection variance to the original is
    # (var'/var)^2 ~= (ab)^2 Imax ~= Imax / 400. So filter any groups where the
    # internal variance is more than 10x what would reasonably be expected after
    # error model correction.
    I = Ih_table.intensities
    mu = Ih_table.Ih_values
    g = Ih_table.inverse_scale_factors
    n_h = Ih_table.sum_in_groups(np.full(Ih_table.size, 1.0))

    group_variances = Ih_table.sum_in_groups(((I / g) - mu) ** 2) / (
        n_h - np.full(n_h.size, 1.0)
    )
    avg_variances = Ih_table.sum_in_groups(Ih_table.variances / (g**2)) / n_h
    ratio = group_variances / avg_variances
    sel = ratio < max(50, (np.max(mu) / 40.0))
    logger.debug(
        f"{sel.size - np.count_nonzero(sel)}/{sel.size} symmetry groups excluded "
        "from error model analysis due to high internal variance"
    )
    Ih_table = Ih_table.select_on_groups(sel)
    n = Ih_table.size
    if n < min_reflections_required:
        raise ValueError(
            "Insufficient reflections (%s < %s) to perform error modelling."
            % (n, min_reflections_required)
        )
    return Ih_table
