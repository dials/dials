"""
Classes that take in a scaler and minimisation parameters and
return the scale factors and derivatives of the scale factors w.r.t.
the parameters
"""

from __future__ import annotations

from scitbx import sparse

from dials.array_family import flex
from dials_scaling_ext import row_multiply


class RefinerCalculator:
    """Class that takes in a scaling_apm and calculates the scale factors
    and derivatives for minimisation."""

    @staticmethod
    def _calc_component_scales_derivatives(apm, block_id):
        """Calculate the scales and derivatives for all components for a given
        block, returning each as a list of values from the components."""
        scales = []
        derivatives = []
        for component in apm.components.values():
            sdc = component["object"].calculate_scales_and_derivatives(block_id)
            scales.append(sdc[0])
            derivatives.append(sdc[1])
        return scales, derivatives

    @staticmethod
    def _calculate_scale_factors(apm, block_id, scales):
        """Calculate the overall scale factor for each reflection from individual
        components."""
        if not scales:
            return apm.constant_g_values[block_id]  # needs to return to set in Ih_table
        multiplied_scale_factors = flex.double(scales[0].size(), 1.0)
        for s in scales:
            multiplied_scale_factors *= s
        if apm.constant_g_values:
            multiplied_scale_factors *= apm.constant_g_values[block_id]
        return multiplied_scale_factors

    @staticmethod
    def _calculate_derivatives(apm, block_id, scales, derivatives_list):
        """Calculate the derivatives matrix."""
        if not scales:
            return sparse.matrix(0, 0)
        if len(scales) == 1:
            # only one active parameter, so don't need to chain rule any derivatives
            # for block_id in range(len(apm.n_obs)):
            return derivatives_list[0]
        derivatives = sparse.matrix(apm.n_obs[block_id], apm.n_active_params)
        col_idx = 0
        for i, d in enumerate(derivatives_list):
            scale_multipliers = flex.double(apm.n_obs[block_id], 1.0)
            for j, s1 in enumerate(scales):
                if i != j:
                    scale_multipliers *= s1
            if apm.constant_g_values:
                scale_multipliers *= apm.constant_g_values[block_id]
            next_deriv = row_multiply(d, scale_multipliers)
            derivatives.assign_block(next_deriv, 0, col_idx)
            col_idx += d.n_cols
        return derivatives

    @classmethod
    def calculate_scales_and_derivatives(cls, apm, block_id):
        """Calculate scale factors and derivatives for minimisation."""
        scales, derivatives = cls._calc_component_scales_derivatives(apm, block_id)
        return (
            cls._calculate_scale_factors(apm, block_id, scales),
            cls._calculate_derivatives(apm, block_id, scales, derivatives),
        )
