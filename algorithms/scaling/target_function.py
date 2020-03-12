"""
This file defines targets for scaling.

These are initialised with a scaler and an active parameter manager,
and have implementations of residual/gradient calculations for
scaling.
"""
from __future__ import absolute_import, division, print_function
from dials.array_family import flex
from dials.algorithms.scaling.scaling_restraints import ScalingRestraintsCalculator
from dials_scaling_ext import row_multiply, calc_dIh_by_dpi, calc_jacobian


class ScalingTarget(object):
    """
    A class to be used by a Scaling Refinery to calculate gradients,
    residuals etc required by the Refinery for minimisation.
    '"""

    _grad_names = ["dI_dp"]
    rmsd_names = ["RMSD_I"]
    rmsd_units = ["a.u"]

    def __init__(self):
        self.rmsd_names = ["RMSD_I"]
        self.rmsd_units = ["a.u"]
        # Quantities to cache each step
        self._rmsds = None
        self.param_restraints = True  # If one tests for restraints and None is
        # returned, then this is set to False and restraints calculations are not
        # attempted for the remainder of the minimisation with this target function.

    def rmsds(self, Ih_table, apm):
        """Calculate RMSDs for the matches. Also calculate R-factors."""
        R = flex.double([])
        n = 0
        for block in Ih_table.blocked_data_list:
            R.extend(flex.pow2(self.calculate_residuals(block)) * block.weights)
            n += block.size
        if self.param_restraints:
            restraints = ScalingRestraintsCalculator.calculate_restraints(apm)
            if restraints:
                R.extend(restraints[0])
            else:
                self.param_restraints = False
        self._rmsds = [(flex.sum(R) / n) ** 0.5]
        return self._rmsds

    @staticmethod
    def achieved():
        """Method required by refinement engine."""
        return False  # implement a method here?

    @staticmethod
    def calculate_residuals(Ih_table):
        """Return the residual vector."""
        R = Ih_table.intensities - (Ih_table.inverse_scale_factors * Ih_table.Ih_values)
        return R

    @staticmethod
    def calculate_gradients(Ih_table):
        """Return a gradient vector on length len(self.apm.x)."""
        gsq = flex.pow2(Ih_table.inverse_scale_factors) * Ih_table.weights
        sumgsq = gsq * Ih_table.h_index_matrix
        prefactor = (
            -2.0
            * Ih_table.weights
            * (
                Ih_table.intensities
                - (Ih_table.Ih_values * Ih_table.inverse_scale_factors)
            )
        )
        dIh = (
            Ih_table.intensities
            - (Ih_table.Ih_values * 2.0 * Ih_table.inverse_scale_factors)
        ) * Ih_table.weights
        dIh_by_dpi = calc_dIh_by_dpi(
            dIh, sumgsq, Ih_table.h_index_matrix, Ih_table.derivatives.transpose()
        )
        term_1 = (prefactor * Ih_table.Ih_values) * Ih_table.derivatives
        term_2 = (
            prefactor * Ih_table.inverse_scale_factors * Ih_table.h_index_matrix
        ) * dIh_by_dpi
        gradient = term_1 + term_2
        return gradient

    @staticmethod
    def calculate_jacobian(Ih_table):
        """Calculate the jacobian matrix, size Ih_table.size by len(self.apm.x)."""
        gsq = flex.pow2(Ih_table.inverse_scale_factors) * Ih_table.weights
        sumgsq = gsq * Ih_table.h_index_matrix
        dIh = (
            Ih_table.intensities
            - (Ih_table.Ih_values * 2.0 * Ih_table.inverse_scale_factors)
        ) * Ih_table.weights
        jacobian = calc_jacobian(
            Ih_table.derivatives.transpose(),
            Ih_table.h_index_matrix,
            Ih_table.Ih_values,
            Ih_table.inverse_scale_factors,
            dIh,
            sumgsq,
        )
        return jacobian

    # The following methods are for adaptlbfgs.
    @classmethod
    def compute_functional_gradients(cls, Ih_table):
        """Return the functional and gradients."""
        resids = cls.calculate_residuals(Ih_table)
        gradients = cls.calculate_gradients(Ih_table)
        weights = Ih_table.weights
        functional = flex.sum(flex.pow2(resids) * weights)
        del Ih_table.derivatives
        return functional, gradients

    def compute_restraints_functional_gradients(self, apm):
        """Return the restrains for functional and gradients."""
        restraints = None
        if self.param_restraints:
            restr = ScalingRestraintsCalculator.calculate_restraints(apm)
            if restr:
                resid_restr = flex.sum(restr[0])  # add to total functional here
                grad_restr = restr[1]
                restraints = [resid_restr, grad_restr]
            else:
                self.param_restraints = False
        return restraints  # list of restraints to add to resid, grads and curvs

    # The following methods are for adaptlstbx (GN/ LM algorithms)
    @classmethod
    def compute_residuals(cls, Ih_table):
        """Return the residuals array and weights."""
        residuals = cls.calculate_residuals(Ih_table)
        weights = Ih_table.weights
        return residuals, weights

    @classmethod
    def compute_residuals_and_gradients(cls, Ih_table):
        """Return the residuals array, jacobian matrix and weights."""
        residuals = cls.calculate_residuals(Ih_table)
        jacobian = cls.calculate_jacobian(Ih_table)
        weights = Ih_table.weights
        Ih_table.derivatives = None
        return residuals, jacobian, weights

    def compute_restraints_residuals_and_gradients(self, apm):
        """Return the restraints for the residuals and jacobian."""
        if self.param_restraints:
            restr = ScalingRestraintsCalculator.calculate_jacobian_restraints(apm)
            if not restr:
                self.param_restraints = False
            return restr
        return None


class ScalingTargetFixedIH(ScalingTarget):
    """An implementation of scaling target for when the scaling is to be
    done against a fixed reference Ih set (i.e scaler is a TargetScaler)
    """

    @staticmethod
    def calculate_gradients(Ih_table):
        rhl = Ih_table.intensities - (
            Ih_table.Ih_values * Ih_table.inverse_scale_factors
        )
        G = -2.0 * rhl * Ih_table.weights * Ih_table.Ih_values * Ih_table.derivatives
        return G

    @staticmethod
    def calculate_jacobian(Ih_table):
        """Calculate the jacobian matrix, size Ih_table.size by len(self.apm.x)."""
        jacobian = row_multiply(Ih_table.derivatives, -1.0 * Ih_table.Ih_values)
        return jacobian
