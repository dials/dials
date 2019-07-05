"""
Error model classes for scaling.
"""
from __future__ import absolute_import, division, print_function
import logging
from math import log, exp
from dials.array_family import flex
from scitbx import sparse
from libtbx.table_utils import simple_table

logger = logging.getLogger("dials.scale")


def get_error_model(error_model_type):
    """Return the correct error model class from a params option."""
    if error_model_type == "basic":
        return BasicErrorModel
    else:
        raise ValueError("Invalid choice of error model: %s" % error_model_type)


class BasicErrorModel(object):
    """
    Object to manage calculation of deviations for an error model.
    """

    min_reflections_required = 250

    def __init__(self, Ih_table, n_bins=10, min_Ih=25.0, min_partiality=0.95):
        logger.info("Initialising an error model for refinement.")
        self.Ih_table = Ih_table
        self.n_bins = n_bins
        self.binning_info = {}
        # First select on initial delta
        self.filter_unsuitable_reflections(
            cutoff=12.0, min_Ih=min_Ih, min_partiality=min_partiality
        )
        self.n_h = self.Ih_table.calc_nh()
        self.sigmaprime = None
        self.delta_hl = None
        self.bin_variances = None
        self._summation_matrix = self.create_summation_matrix()
        self._bin_counts = flex.double(self.Ih_table.size, 1.0) * self.summation_matrix
        self.weights = self._bin_counts ** 0.5
        self.refined_parameters = [1.0, 0.0]

    def __str__(self):
        a = abs(self.refined_parameters[0])
        b = abs(self.refined_parameters[1])
        ISa = "%.3f" % (1.0 / (b * a)) if (b * a) > 0 else "Unable to estimate"
        return "\n".join(
            (
                "",
                "Error model details:",
                "  Type: basic",
                "  Current parameters: a = %.5f, b = %.5f" % (a, b),
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
        """Output a summary of model minimisation to the logger."""
        header = ["Intensity range (<Ih>)", "n_refl", "variance(norm_dev)"]
        rows = []
        bin_bounds = ["%.2f" % i for i in self.binning_info["bin_boundaries"]]
        for i, (bin_var, n_refl) in enumerate(
            zip(self.binning_info["bin_variances"], self.binning_info["refl_per_bin"])
        ):
            rows.append(
                [
                    bin_bounds[i] + " - " + bin_bounds[i + 1],
                    str(n_refl),
                    str(round(bin_var, 3)),
                ]
            )
        st = simple_table(rows, header)
        logger.info(
            "\n".join(
                (
                    "Intensity bins used during error model refinement:",
                    st.format(),
                    "variance(norm_dev) expected to be ~ 1 for each bin.",
                    "",
                )
            )
        )

    @property
    def summation_matrix(self):
        """A sparse matrix to allow summation over intensity groups."""
        return self._summation_matrix

    @property
    def bin_counts(self):
        """An array of the number of intensities assigned to each bin."""
        return self._bin_counts

    def filter_unsuitable_reflections(self, cutoff, min_Ih, min_partiality):
        """Do a first pass to calculate delta_hl and filter out the largest
        deviants, so that the error model is not misled by these and instead
        operates on the central ~90% of the data. Also choose reflection groups
        with n_h > 1, as these have deltas of zero by definition and will bias
        the variance calculations. Also, only use groups where <Ih> > 2.0, as
        the assumptions of normally distributed deltas will not hold for low
        <Ih>."""
        self.n_h = self.Ih_table.calc_nh()
        self.sigmaprime = self.calc_sigmaprime([1.0, 0.0])
        delta_hl = self.calc_deltahl()
        # make sure the fit isn't misled by extreme values
        sel = flex.abs(delta_hl) < cutoff
        if "partiality" in self.Ih_table.Ih_table:
            sel &= self.Ih_table.Ih_table["partiality"] > min_partiality
        self.Ih_table = self.Ih_table.select(sel)

        n = self.Ih_table.size
        sum_I_over_var = (
            self.Ih_table.intensities / self.Ih_table.variances
        ) * self.Ih_table.h_index_matrix
        n_per_group = flex.double(n, 1) * self.Ih_table.h_index_matrix
        avg_I_over_var = sum_I_over_var / n_per_group
        sel = avg_I_over_var > 0.85
        self.Ih_table = self.Ih_table.select_on_groups(sel)
        self.n_h = self.Ih_table.calc_nh()
        scaled_Ih = self.Ih_table.Ih_values * self.Ih_table.inverse_scale_factors
        # need a scaled min_Ih, where can reasonably expect norm distribution
        # (use min_Ih=25 by default, sigma ~ 5)
        sel2 = scaled_Ih > min_Ih
        # can't calculate a true deviation for groups of 1
        sel3 = self.n_h > 1.0
        sel4 = self.Ih_table.intensities > 0.001
        # don't want to include weaker reflections where the background adds
        # significantly to the variances, as these would no longer be normally
        # distributed and skew the fit.
        self.Ih_table = self.Ih_table.select(sel2 & sel3 & sel4)
        n = self.Ih_table.size
        if n < self.min_reflections_required:
            raise ValueError(
                "Insufficient reflections (%s) to perform error modelling." % n
            )
        self.n_h = self.Ih_table.calc_nh()
        # now make sure any left also have n > 1
        sel = self.n_h > 1.0
        self.Ih_table = self.Ih_table.select(sel)
        self.n_h = self.Ih_table.calc_nh()

    def calc_sigmaprime(self, x):
        """Calculate the error from the model."""
        sigmaprime = (
            x[0]
            * ((self.Ih_table.variances) + ((x[1] * self.Ih_table.intensities) ** 2))
            ** 0.5
        ) / self.Ih_table.inverse_scale_factors
        return sigmaprime

    def calc_deltahl(self):
        """Calculate the normalised deviations from the model."""
        I_hl = self.Ih_table.intensities
        g_hl = self.Ih_table.inverse_scale_factors
        I_h = self.Ih_table.Ih_values
        prefactor = ((self.n_h - flex.double(self.n_h.size(), 1.0)) / self.n_h) ** 0.5
        delta_hl = prefactor * ((I_hl / g_hl) - I_h) / self.sigmaprime
        return delta_hl

    def update_for_minimisation(self, x):
        """"Calculate the updated quantites."""
        self.sigmaprime = self.calc_sigmaprime(x)
        self.delta_hl = self.calc_deltahl()
        self.bin_variances = self.calculate_bin_variances()

    def create_summation_matrix(self):
        """"Create a summation matrix to allow sums into intensity bins.

        This routine attempts to bin into bins equally spaced in log(intensity),
        to give a representative sample across all intensities. To avoid
        undersampling, it is required that there are at least 100 reflections
        per intensity bin unless there are very few reflections."""
        n = self.Ih_table.size
        if n < self.min_reflections_required:
            raise ValueError(
                "Insufficient reflections (%s) to perform error modelling." % n
            )
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
        self.binning_info["refl_per_bin"] = []

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
        new_sum_matrix = sparse.matrix(summation_matrix.n_rows, n_new_cols)
        next_col = 0
        for i, col in enumerate(summation_matrix.cols()):
            if i not in cols_to_del:
                new_sum_matrix[:, next_col] = col
                next_col += 1
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

    def update_variances(self, variances, intensities):
        """Use the error model parameter to calculate new values for the variances."""
        new_variance = (self.refined_parameters[0] ** 2) * (
            variances + ((self.refined_parameters[1] * intensities) ** 2)
        )
        return new_variance

    def clear_Ih_table(self):
        """Delete the Ih_table, to free memory."""
        del self.Ih_table
