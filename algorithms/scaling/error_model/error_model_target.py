"""
Definition of the target function for error model minimisation.
"""

from __future__ import absolute_import, division, print_function
from dials.array_family import flex


class ErrorModelTarget(object):

    """Error model target for finding slope of norm distribution.
    (i.e. the 'a' parameter of the basic errror model)"""

    _grad_names = ["d(deltahl)_dp"]
    rmsd_names = ["RMSD_deltahl"]
    rmsd_units = ["a.u"]

    def __init__(self, error_model, parameterisation):
        self.error_model = error_model
        self.parameterisation = parameterisation
        # intercept + gradient
        # Quantities to cache each step
        self._rmsds = None

    def predict(self):
        """Do the next step of the prediction."""
        self.error_model.update_for_minimisation(self.parameterisation)

    def get_num_matches(self):
        """Get the number of reflections."""
        return self.error_model.n_refl

    def rmsds(self):
        """calculate unweighted RMSDs for the matches"""
        # cache rmsd calculation for achieved test
        R = self.calculate_residuals()
        n = R.size()
        R.extend(flex.double([self.compute_restraints_functional_gradients()[0]]))
        self._rmsds = [(flex.sum(R) / n) ** 0.5]
        return self._rmsds

    def achieved(self):
        """Method required for refinement engine."""
        return False  # implement a method here?

    # The following methods are for adaptlbfgs.
    def compute_functional_gradients(self):
        """Compute the functional and gradients vector."""
        return flex.sum(self.calculate_residuals()), self.calculate_gradients()

    def compute_restraints_functional_gradients(self):
        """Compute the restraints for the functional and gradients."""
        return [0.0, 0.0]


class ErrorModelTargetA(ErrorModelTarget):

    """Target to minimise the 'a' component of the basic error model."""

    def __init__(self, error_model, parameterisation):
        super(ErrorModelTargetA, self).__init__(
            error_model.components["a"], parameterisation
        )
        self.parameterisation.set_active_parameter("a")

    def predict(self):
        """Do the next step of the prediction."""
        pass

    def calculate_residuals(self):
        """Return the residual vector"""
        x = self.parameterisation.x
        R = (self.error_model.sortedy - (x[1] * self.error_model.sortedx) - x[0]) ** 2
        return R

    def calculate_gradients(self):
        "calculate the gradient vector"
        x = self.parameterisation.x
        R = self.error_model.sortedy - (x[1] * self.error_model.sortedx) - x[0]
        gradient = flex.double(
            [-2.0 * flex.sum(R), -2.0 * flex.sum(R * self.error_model.sortedx)]
        )
        return gradient


class ErrorModelTargetB(ErrorModelTarget):

    """Target to minimise the 'b' component of the basic error model."""

    def __init__(self, error_model, parameterisation):
        super(ErrorModelTargetB, self).__init__(
            error_model.components["b"], parameterisation
        )
        self.parameterisation.set_active_parameter("b")

    def calculate_residuals(self):
        """Return the residual vector"""
        bin_vars = self.error_model.bin_variances
        R = (
            (
                (flex.double(bin_vars.size(), 0.5) - bin_vars) ** 2
                + (1.0 / bin_vars)
                - flex.double(bin_vars.size(), 1.25)
            )
            * self.error_model.weights
            / flex.sum(self.error_model.weights)
        )
        return R

    def calculate_gradients(self):
        "calculate the gradient vector"
        a, b = self.parameterisation.a, self.parameterisation.x[0]
        Ih_table = self.error_model.Ih_table
        I_hl = Ih_table.intensities
        g_hl = Ih_table.inverse_scale_factors
        weights = self.error_model.weights
        bin_vars = self.error_model.bin_variances
        sum_matrix = self.error_model.summation_matrix
        bin_counts = self.error_model.bin_counts
        dsig_dc = (
            b * (I_hl ** 2) * (a ** 2) / (self.error_model.sigmaprime * (g_hl ** 2))
        )
        ddelta_dsigma = -1.0 * self.error_model.delta_hl / self.error_model.sigmaprime
        deriv = ddelta_dsigma * dsig_dc
        dphi_by_dvar = -2.0 * (
            flex.double(bin_vars.size(), 0.5)
            - bin_vars
            + (1.0 / (2.0 * (bin_vars ** 2)))
        )
        term1 = 2.0 * self.error_model.delta_hl * deriv * sum_matrix
        term2a = self.error_model.delta_hl * sum_matrix
        term2b = deriv * sum_matrix
        grad = dphi_by_dvar * (
            (term1 / bin_counts) - (2.0 * term2a * term2b / (bin_counts ** 2))
        )
        gradients = flex.double([flex.sum(grad * weights) / flex.sum(weights)])
        return gradients
