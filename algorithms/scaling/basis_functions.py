"""
Classes that take in a scaler and minimisation parameters and
return the scale factors and derivatives of the scale factors w.r.t.
the parameters
"""
from __future__ import absolute_import, division, print_function

from dials.array_family import flex
from dials_scaling_ext import row_multiply
from scitbx import sparse


class basis_function(object):
    """Class that takes in a scaling_apm and calcuates the scale factors,
    derivatives and optionally curvatures for minimisation."""

    def calc_component_scales_derivatives(self, apm, block_id):
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
    def calculate_scale_factors(apm, block_id, scales):
        """Calculate the overall scale factor for each reflection from individual
        components."""
        if not scales:
            return None
        multiplied_scale_factors = flex.double(scales[0].size(), 1.0)
        for s in scales:
            multiplied_scale_factors *= s
        if apm.constant_g_values:
            multiplied_scale_factors *= apm.constant_g_values[block_id]
        return multiplied_scale_factors

    @staticmethod
    def calculate_derivatives(apm, block_id, scales, derivatives_list):
        """Calculate the derivatives matrix."""
        if not scales:
            return None
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

    def calculate_scales_and_derivatives(self, apm, block_id):
        """Calculate and return scale factors, derivatives and optionally
        curvatures to be used in minimisation."""
        scales, derivatives = self.calc_component_scales_derivatives(apm, block_id)
        return (
            self.calculate_scale_factors(apm, block_id, scales),
            self.calculate_derivatives(apm, block_id, scales, derivatives),
        )
