"""
This file defines targets for scaling.

These are initialised with a scaler and an active parameter manager,
and have implementations of residual/gradient calculations for
scaling.
"""

from __future__ import annotations

import numpy as np

from dxtbx import flumpy

from dials.algorithms.scaling.scaling_restraints import ScalingRestraintsCalculator
from dials.array_family import flex
from dials_scaling_ext import calc_dIh_by_dpi, calc_jacobian, row_multiply


class ScalingTarget:
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
        R = 0
        n = 0
        for block in Ih_table.blocked_data_list:
            R += np.sum(np.square(self.calculate_residuals(block)) * block.weights)
            n += block.size
        if self.param_restraints:
            restraints = ScalingRestraintsCalculator.calculate_restraints(apm)
            if restraints:
                R += np.sum(restraints[0])
            else:
                self.param_restraints = False
        self._rmsds = [(R / n) ** 0.5]
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
        gsq = np.square(Ih_table.inverse_scale_factors) * Ih_table.weights
        sumgsq = flumpy.from_numpy(Ih_table.sum_in_groups(gsq))
        prefactor = (
            -2.0
            * Ih_table.weights
            * (
                Ih_table.intensities
                - (Ih_table.Ih_values * Ih_table.inverse_scale_factors)
            )
        )
        dIh = flumpy.from_numpy(
            (
                Ih_table.intensities
                - (Ih_table.Ih_values * 2.0 * Ih_table.inverse_scale_factors)
            )
            * Ih_table.weights
        )
        dIh_by_dpi = calc_dIh_by_dpi(
            dIh, sumgsq, Ih_table.h_index_matrix, Ih_table.derivatives.transpose()
        )
        term_1 = (
            flumpy.from_numpy(prefactor * Ih_table.Ih_values) * Ih_table.derivatives
        )
        term_2 = (
            flumpy.from_numpy(
                Ih_table.sum_in_groups(prefactor * Ih_table.inverse_scale_factors)
            )
            * dIh_by_dpi
        )
        gradient = term_1 + term_2
        return gradient

    @staticmethod
    def calculate_jacobian(Ih_table):
        """Calculate the jacobian matrix, size Ih_table.size by len(self.apm.x)."""
        gsq = np.square(Ih_table.inverse_scale_factors) * Ih_table.weights
        sumgsq = Ih_table.sum_in_groups(gsq)
        dIh = (
            Ih_table.intensities
            - (Ih_table.Ih_values * 2.0 * Ih_table.inverse_scale_factors)
        ) * Ih_table.weights
        jacobian = calc_jacobian(
            Ih_table.derivatives.transpose(),
            Ih_table.h_index_matrix,
            flumpy.from_numpy(Ih_table.Ih_values),
            flumpy.from_numpy(Ih_table.inverse_scale_factors),
            flumpy.from_numpy(dIh),
            flumpy.from_numpy(sumgsq),
        )
        return jacobian

    # The following methods are for adaptlbfgs.
    @classmethod
    def compute_functional_gradients(cls, Ih_table):
        """Return the functional and gradients."""
        resids = cls.calculate_residuals(Ih_table)
        gradients = cls.calculate_gradients(Ih_table)
        weights = Ih_table.weights
        functional = np.sum(np.square(resids) * weights)
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
        return flumpy.from_numpy(residuals), flumpy.from_numpy(weights)

    @classmethod
    def compute_residuals_and_gradients(cls, Ih_table):
        """Return the residuals array, jacobian matrix and weights."""
        residuals = cls.calculate_residuals(Ih_table)
        jacobian = cls.calculate_jacobian(Ih_table)
        weights = Ih_table.weights
        Ih_table.derivatives = None
        return flumpy.from_numpy(residuals), jacobian, flumpy.from_numpy(weights)

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
        G = (
            flumpy.from_numpy(-2.0 * rhl * Ih_table.weights * Ih_table.Ih_values)
            * Ih_table.derivatives
        )
        return G

    @staticmethod
    def calculate_jacobian(Ih_table):
        """Calculate the jacobian matrix, size Ih_table.size by len(self.apm.x)."""
        jacobian = row_multiply(
            Ih_table.derivatives, flumpy.from_numpy(-1.0 * Ih_table.Ih_values)
        )
        return jacobian
