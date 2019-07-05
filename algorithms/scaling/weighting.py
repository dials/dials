"""
Weighting scheme definitions for scaling.
"""
from __future__ import absolute_import, division, print_function
from dials.array_family import flex


def get_weighting_scheme(Ih_table, weighting_scheme):
    """Return the appropriate weighting scheme from the params option."""
    if weighting_scheme == "invvar" or weighting_scheme is None:
        return WeightingBase(Ih_table)
    if weighting_scheme == "unity":
        return UnityWeights(Ih_table)
    if weighting_scheme == "GM":
        return GemanMcClureWeights(Ih_table)
    if weighting_scheme == "cauchy":
        return CauchyWeights(Ih_table)
    if weighting_scheme == "huber":
        return HuberWeights(Ih_table)
    else:
        raise ValueError("Invalid choice of weighting scheme: %s" % weighting_scheme)


class WeightingBase(object):
    """Base class that defines the properties of a Scaling Weights object."""

    weighting_scheme = "invvar"

    def __init__(self, Ih_table):
        self.Ih_table = Ih_table

    def calculate_initial_weights(self):
        """Set the initial weight, to allow first calculation of Ih."""
        self.Ih_table.weights = 1.0 / self.Ih_table.variances

    def apply_iterative_weights(self):
        """Update the weights, if relevant."""
        pass


class UnityWeights(WeightingBase):
    """Class for unity weights."""

    weighting_scheme = "unity"

    def calculate_initial_weights(self):
        self.Ih_table.weights = flex.double(self.Ih_table.size, 1.0)


class GemanMcClureWeights(WeightingBase):
    """Class for Geman-McClure weighting scheme (M-estimation)."""

    weighting_scheme = "GM"

    def apply_iterative_weights(self):
        e_i = (
            self.Ih_table.intensities
            - (self.Ih_table.inverse_scale_factors * self.Ih_table.Ih_values)
        ) / (self.Ih_table.inverse_scale_factors * self.Ih_table.Ih_values)
        self.Ih_table.weights = 1.0 / ((1.0 + (e_i ** 2)) ** 2)


class CauchyWeights(WeightingBase):
    """Class for Geman-McClure weighting scheme (M-estimation)."""

    weighting_scheme = "cauchy"

    c = 2.385  # Value for 95% efficiency

    def apply_iterative_weights(self):
        e_i = (
            self.Ih_table.intensities
            - (self.Ih_table.inverse_scale_factors * self.Ih_table.Ih_values)
        ) / (self.Ih_table.inverse_scale_factors * self.Ih_table.Ih_values)
        self.Ih_table.weights = 1.0 / (1.0 + ((e_i / self.c) ** 2))


class HuberWeights(WeightingBase):
    """Class for Geman-McClure weighting scheme (M-estimation)."""

    weighting_scheme = "huber"

    c = 1.345  # Value for 95% efficiency

    def apply_iterative_weights(self):
        e_i = (
            self.Ih_table.intensities
            - (self.Ih_table.inverse_scale_factors * self.Ih_table.Ih_values)
        ) / (self.Ih_table.inverse_scale_factors * self.Ih_table.Ih_values)
        abs_ei = (e_i ** 2) ** 0.5
        self.Ih_table.weights = flex.double(self.Ih_table.size, 1.0)
        sel = abs_ei > self.c
        sel_abs_ei = abs_ei.select(sel)
        sel_weight_fn = self.c / sel_abs_ei
        self.Ih_table.weights.set_selected(sel, sel_weight_fn)
