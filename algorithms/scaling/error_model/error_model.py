"""
Error model classes for scaling.
"""
from __future__ import absolute_import, division, print_function
import logging
from math import log, exp
from dials.util import tabulate
from dials.array_family import flex
from scitbx import sparse
from scitbx.math.distributions import normal_distribution

logger = logging.getLogger("dials")


def get_error_model(error_model_type):
    """Return the correct error model class from a params option."""
    if error_model_type == "basic":
        return BasicErrorModel
    else:
        raise ValueError("Invalid choice of error model: %s" % error_model_type)


def calc_sigmaprime(x, Ih_table):
    """Calculate the error from the model."""
    sigmaprime = (
        x[0] * ((Ih_table.variances) + ((x[1] * Ih_table.intensities) ** 2)) ** 0.5
    ) / Ih_table.inverse_scale_factors
    return sigmaprime


def calc_deltahl(Ih_table, n_h, sigmaprime):
    """Calculate the normalised deviations from the model."""
    I_hl = Ih_table.intensities
    g_hl = Ih_table.inverse_scale_factors
    I_h = Ih_table.Ih_values
    prefactor = ((n_h - flex.double(n_h.size(), 1.0)) / n_h) ** 0.5
    delta_hl = prefactor * ((I_hl / g_hl) - I_h) / sigmaprime
    return delta_hl


class BasicErrorModelParameterisation(object):

    """Parameter manager for error model minimisation."""

    def __init__(self, basic_error_model):
        p = basic_error_model.parameters
        self.a = p[0]
        self.b = p[1]
        self.x = None
        self.active_parameterisation = None
        if "a" in basic_error_model.free_components:
            self.set_active_parameter("a")
        elif "b" in basic_error_model.free_components:
            self.set_active_parameter("b")

    def set_param_vals(self, x):
        """method for refinement engine access"""
        self.x = x
        if self.active_parameterisation == "a":
            self.a = self.x[1]
        elif self.active_parameterisation == "b":
            self.b = self.x[0]

    def get_param_vals(self):
        """method for refinement engine access"""
        return self.x

    def set_active_parameter(self, name):
        """Set the active parameter being minimised."""
        if name == "b":
            self.active_parameterisation = "b"
            self.x = [self.b]
        elif name == "a":
            self.active_parameterisation = "a"
            self.x = [0.0, self.a]


class BasicErrorModel(object):

    """Two component error model."""

    min_reflections_required = 250

    def __init__(self, Ih_table, params):

        self.components = {}
        self.free_components = []
        self.filtered_Ih_table = None
        a, b = (1.0, 0.02)

        error_params = params.weighting.error_model
        if error_params.basic.a:
            a = error_params.basic.a
        if error_params.basic.b:
            b = error_params.basic.b

        self.components = {
            "a": BasicErrorModelA(a),
            "b": BasicErrorModelB(b, error_params.n_bins),
        }
        try:
            self.filtered_Ih_table = self.filter_unsuitable_reflections(
                Ih_table,
                cutoff=12.0,
                min_Ih=error_params.min_Ih,
                min_partiality=params.reflection_selection.min_partiality,
            )
        except ValueError as e:
            """No free components set, no minimisation."""
            logger.info(e)
        else:
            # See if suitable data for minimisation.
            try:
                self.components["a"].add_data(self.filtered_Ih_table)
            except ValueError:
                pass
            else:
                if params:
                    if not error_params.basic.a:
                        self.free_components.append("a")
            self.components["b"].add_data(self.filtered_Ih_table)
            if params:
                if not error_params.basic.b:
                    self.free_components.append("b")

    @property
    def parameters(self):
        """A list of the model parameters."""
        return [self.components["a"].parameters[0], self.components["b"].parameters[0]]

    @parameters.setter
    def parameters(self, parameters):
        assert len(parameters) == 2
        self.components["a"].parameters = [parameters[0]]
        self.components["b"].parameters = [parameters[1]]

    @classmethod
    def filter_unsuitable_reflections(cls, Ih_table, cutoff, min_Ih, min_partiality):
        """Do a first pass to calculate delta_hl and filter out the largest
        deviants, so that the error model is not misled by these and instead
        operates on the central ~90% of the data. Also choose reflection groups
        with n_h > 1, as these have deltas of zero by definition and will bias
        the variance calculations. Also, only use groups where <Ih> > 25.0, as
        the assumptions of normally distributed deltas will not hold for low
        <Ih>."""
        n_h = Ih_table.calc_nh()
        sigmaprime = calc_sigmaprime([1.0, 0.0], Ih_table)
        delta_hl = calc_deltahl(Ih_table, n_h, sigmaprime)
        # make sure the fit isn't misled by extreme values
        sel = flex.abs(delta_hl) < cutoff
        if "partiality" in Ih_table.Ih_table:
            sel &= Ih_table.Ih_table["partiality"] > min_partiality
        Ih_table = Ih_table.select(sel)

        n = Ih_table.size
        sum_I_over_var = (
            Ih_table.intensities / Ih_table.variances
        ) * Ih_table.h_index_matrix
        n_per_group = flex.double(n, 1) * Ih_table.h_index_matrix
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
        if n < cls.min_reflections_required:
            raise ValueError(
                "Insufficient reflections (%s < %s) to perform error modelling."
                % (n, cls.min_reflections_required)
            )
        n_h = Ih_table.calc_nh()
        # now make sure any left also have n > 1
        sel = n_h > 1.0
        Ih_table = Ih_table.select(sel)
        return Ih_table

    def update_parameters(self, parameterisation):
        """Update the state of the error model."""
        self.components["b"].parameters = [parameterisation.b]
        if self.components["b"].Ih_table:
            self.components["b"].update_for_minimisation(parameterisation)
        self.components["a"].parameters = [parameterisation.a]
        if "a" in self.free_components:
            self.components["a"].reinitialise([parameterisation.a, parameterisation.b])

    def update_parameters_after_minimisation(self, parameterisation):
        """Update the state of the error model and prepare for next minimisation."""
        if parameterisation.active_parameterisation == "a":
            self.components["a"].parameters = [parameterisation.a]
            # just refined a, so need to prepare b.
            if self.components["b"].Ih_table:
                self.components["b"].update_for_minimisation(parameterisation)
        elif parameterisation.active_parameterisation == "b":
            # update n
            self.components["b"].parameters = [parameterisation.b]
            # just refined b, so prepare a - need to reinitialise with 1.0,
            # so that can determine new slope.
            if "a" in self.free_components:
                self.components["a"].reinitialise([1.0, parameterisation.b])

    def __str__(self):
        a = abs(self.parameters[0])
        b = abs(self.parameters[1])
        ISa = "%.3f" % (1.0 / (b * a)) if (b * a) > 0 else "Unable to estimate"
        return "\n".join(
            (
                "",
                "Error model details:",
                "  Type: basic",
                "  Parameters: a = %.5f, b = %.5f" % (a, b),
                "  Error model formula: "
                + u"\u03C3"
                + "'"
                + u"\xb2"
                + " = a"
                + u"\xb2"
                + "("
                + u"\u03C3\xb2"
                " + (bI)" + u"\xb2" + ")",
                "  estimated I/sigma asymptotic limit: %s" % ISa,
                "",
            )
        )

    def minimisation_summary(self):
        """Print a summary of minimisation."""
        if self.components["b"].Ih_table:
            return self.components["b"].minimisation_summary()

    def update_variances(self, variances, intensities):
        """Use the error model parameter to calculate new values for the variances."""
        new_variance = (self.parameters[0] ** 2) * (
            variances + ((self.parameters[1] * intensities) ** 2)
        )
        return new_variance

    def clear_Ih_table(self):
        """Delete the Ih_table, to free memory."""
        del self.components["a"].Ih_table
        del self.components["b"].Ih_table
        del self.filtered_Ih_table


class BasicErrorModelB(object):
    """
    Object to manage calculation of deviations for an error model.
    """

    min_reflections_required = 250

    def __init__(self, b=0.0, n_bins=10):
        self.parameters = [b]
        self.Ih_table = None
        self.n_h = None
        self.sigmaprime = None
        self.bin_variances = None
        self._summation_matrix = None
        self.weights = None
        self.delta_hl = None
        self.n_bins = n_bins
        self.binning_info = {
            "initial_variances": [],
            "bin_boundaries": [],
            "mean_intensities": [],
            "bin_variances": [],
            "refl_per_bin": [],
            "n_reflections": None,
        }

    def add_data(self, Ih_table):
        """Try to add data to the component."""
        self.Ih_table = Ih_table
        # First select on initial delta
        self.n_h = self.Ih_table.calc_nh()
        self.sigmaprime = calc_sigmaprime([1.0, 0.0], self.Ih_table)
        self._summation_matrix = self.create_summation_matrix()
        self.weights = flex.double(self.binning_info["mean_intensities"])
        self.delta_hl = calc_deltahl(self.Ih_table, self.n_h, self.sigmaprime)
        self.bin_variances = self.calculate_bin_variances()
        self.binning_info["initial_variances"] = self.binning_info["bin_variances"]

    def minimisation_summary(self):
        """Generate a summary of the model minimisation for output."""
        header = [
            "Intensity range (<Ih>)",
            "n_refl",
            "Uncorrected variance",
            "Corrected variance",
        ]
        rows = []
        bin_bounds = ["%.2f" % i for i in self.binning_info["bin_boundaries"]]
        for i, (initial_var, bin_var, n_refl) in enumerate(
            zip(
                self.binning_info["initial_variances"],
                self.binning_info["bin_variances"],
                self.binning_info["refl_per_bin"],
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

    @property
    def summation_matrix(self):
        """A sparse matrix to allow summation over intensity groups."""
        return self._summation_matrix

    @property
    def bin_counts(self):
        """An array of the number of intensities assigned to each bin."""
        return self.binning_info["refl_per_bin"]

    @property
    def n_refl(self):
        """The number of reflections being used in minimisation"""
        return self.Ih_table.size

    def update_for_minimisation(self, parameterisation):
        """"Calculate the updated quantites."""
        self.sigmaprime = calc_sigmaprime(
            [parameterisation.a, parameterisation.b], self.Ih_table
        )
        self.delta_hl = calc_deltahl(self.Ih_table, self.n_h, self.sigmaprime)
        self.bin_variances = self.calculate_bin_variances()

    def create_summation_matrix(self):
        """"Create a summation matrix to allow sums into intensity bins.

        This routine attempts to bin into bins equally spaced in log(intensity),
        to give a representative sample across all intensities. To avoid
        undersampling, it is required that there are at least 100 reflections
        per intensity bin unless there are very few reflections."""
        n = self.Ih_table.size
        self.binning_info["n_reflections"] = n
        summation_matrix = sparse.matrix(n, self.n_bins)
        Ih = self.Ih_table.Ih_values * self.Ih_table.inverse_scale_factors
        size_order = flex.sort_permutation(Ih, reverse=True)
        Imax = max(Ih)
        Imin = max(1.0, min(Ih))  # avoid log issues
        spacing = (log(Imax) - log(Imin)) / float(self.n_bins)
        boundaries = [Imax] + [
            exp(log(Imax) - (i * spacing)) for i in range(1, self.n_bins + 1)
        ]
        boundaries[-1] = min(Ih) - 0.01
        self.binning_info["bin_boundaries"] = boundaries
        self.binning_info["refl_per_bin"] = flex.double()

        n_cumul = 0
        if Ih.size() > 100 * self.min_reflections_required:
            self.min_reflections_required = int(Ih.size() / 100.0)
        min_per_bin = min(self.min_reflections_required, int(n / (3.0 * self.n_bins)))
        for i in range(len(boundaries) - 1):
            maximum = boundaries[i]
            minimum = boundaries[i + 1]
            sel1 = Ih <= maximum
            sel2 = Ih > minimum
            sel = sel1 & sel2
            isel = sel.iselection()
            n_in_bin = isel.size()
            if n_in_bin < min_per_bin:  # need more in this bin
                m = n_cumul + min_per_bin
                if m < n:  # still some refl left to use
                    idx = size_order[m]
                    intensity = Ih[idx]
                    boundaries[i + 1] = intensity
                    minimum = boundaries[i + 1]
                    sel = sel1 & (Ih > minimum)
                    isel = sel.iselection()
                    n_in_bin = isel.size()
            self.binning_info["refl_per_bin"].append(n_in_bin)
            for j in isel:
                summation_matrix[j, i] = 1
            n_cumul += n_in_bin
        cols_to_del = []
        for i, col in enumerate(summation_matrix.cols()):
            if col.non_zeroes < min_per_bin - 5:
                cols_to_del.append(i)
        n_new_cols = summation_matrix.n_cols - len(cols_to_del)
        if n_new_cols == self.n_bins:
            for i in range(len(boundaries) - 1):
                maximum = boundaries[i]
                minimum = boundaries[i + 1]
                sel1 = Ih <= maximum
                sel2 = Ih > minimum
                sel = sel1 & sel2
                m = flex.mean(Ih.select(sel))
                self.binning_info["mean_intensities"].append(m)
            return summation_matrix
        new_sum_matrix = sparse.matrix(summation_matrix.n_rows, n_new_cols)
        next_col = 0
        refl_per_bin = flex.double()
        new_bounds = []
        for i, col in enumerate(summation_matrix.cols()):
            if i not in cols_to_del:
                new_sum_matrix[:, next_col] = col
                next_col += 1
                new_bounds.append(boundaries[i])
                refl_per_bin.append(self.binning_info["refl_per_bin"][i])
        self.binning_info["refl_per_bin"] = refl_per_bin
        new_bounds.append(boundaries[-1])
        self.binning_info["bin_boundaries"] = new_bounds
        for i in range(len(new_bounds) - 1):
            maximum = new_bounds[i]
            minimum = new_bounds[i + 1]
            sel1 = Ih <= maximum
            sel2 = Ih > minimum
            sel = sel1 & sel2
            m = flex.mean(Ih.select(sel))
            self.binning_info["mean_intensities"].append(m)
        return new_sum_matrix

    def calculate_bin_variances(self):
        """Calculate the variance of each bin."""
        sum_deltasq = (self.delta_hl ** 2) * self.summation_matrix
        sum_delta_sq = (self.delta_hl * self.summation_matrix) ** 2
        bin_vars = (sum_deltasq / self.bin_counts) - (
            sum_delta_sq / (self.bin_counts ** 2)
        )
        self.binning_info["bin_variances"] = bin_vars
        return bin_vars


class BasicErrorModelA(object):

    """The pieces needed for determining the 'a' error model component."""

    def __init__(self, a=1.0):
        self.parameters = [a]
        self.Ih_table = None
        self.sortedx = None
        self.sortedy = None

    def add_data(self, Ih_table):
        """Try to add data to component."""
        sel = (Ih_table.intensities / (Ih_table.variances ** 0.5)) < 20.0
        self.Ih_table = Ih_table.select(sel)
        if not self.Ih_table.size:
            raise ValueError(
                "No suitable reflections remain for error model 'a' parameter minimisation."
            )
        self.calc_sortedxy([1.0, 0.02])  # start with a = 1.0, b= 0.02

    @property
    def n_refl(self):
        """The number of reflections being used for minimisation."""
        return self.Ih_table.size

    def reinitialise(self, params):
        """Calculate the normalised deltas."""
        self.calc_sortedxy(params)

    def calc_sortedxy(self, params):
        """Sort the x,y data."""
        sigmaprime = calc_sigmaprime(params, self.Ih_table)
        delta_hl = calc_deltahl(self.Ih_table, self.Ih_table.calc_nh(), sigmaprime)
        norm = normal_distribution()
        n = len(delta_hl)
        if n <= 10:
            a = 3 / 8
        else:
            a = 0.5
        self.sortedy = flex.sorted(flex.double(delta_hl))
        self.sortedx = flex.double(
            [norm.quantile((i + 1 - a) / (n + 1 - (2 * a))) for i in range(n)]
        )
        central_sel = (self.sortedx < 1.5) & (self.sortedx > -1.5)
        self.sortedy = self.sortedy.select(central_sel)
        self.sortedx = self.sortedx.select(central_sel)
