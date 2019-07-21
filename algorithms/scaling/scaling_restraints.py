"""
Restraints manager classes for scaling.
"""
from __future__ import absolute_import, division, print_function
from scitbx import sparse
from dials.array_family import flex


class MultiScalingRestraints(object):
    """A class to calculate restraints for scaling for one or more datasets,
    by composition of the restraints of the individual datasets. The methods
    require a multi active_parameter_manager."""

    def __init__(self):
        self.single_restraints_calculator = ScalingRestraints()

    def calculate_restraints(self, apm):
        """Calculate restraints for multi-dataset scaling, using a
        multi active_parameter_manager. Return None if no restraints, else return
        a restraints vector of length n_restrained_parameters and a gradient vector
        of length n_total_parameters."""
        R = flex.double([])
        G = flex.double([])
        for single_apm in apm.apm_list:
            restr = self.single_restraints_calculator.calculate_restraints(single_apm)
            if restr:
                R.extend(restr[0])
                G.extend(restr[1])
            else:
                G.extend(flex.double(single_apm.n_active_params, 0.0))
        if R:
            return [R, G]
        return None

    def calculate_jacobian_restraints(self, apm):
        """Calculate jacobian restraints for multi-dataset scaling, using a
        multi active_parameter_manager. First, the restraints are calculated - if
        not None is returned, then one restraints vector, jacobian matrix and
        weights vector is composed for the multiple datasets, else None is returned.
        The jacobian restraints matrix is of size n_restrained_parameters x
        n_parameters (across all datasets), while the residuals and weights vector
        are of length n_restrainted_parameters."""
        residual_restraints = self.calculate_restraints(apm)
        if residual_restraints:
            n_restraints = residual_restraints[0].size()
            weights = flex.double([])
            restraints_vector = flex.double([])
            jacobian = sparse.matrix(n_restraints, apm.n_active_params)
            cumul_restr_pos = 0
            for i, single_apm in enumerate(apm.apm_list):
                restraints = self.single_restraints_calculator.calculate_jacobian_restraints(
                    single_apm
                )
                if restraints:
                    jacobian.assign_block(
                        restraints[1], cumul_restr_pos, apm.apm_data[i]["start_idx"]
                    )
                    cumul_restr_pos += restraints[1].n_rows
                    restraints_vector.extend(restraints[0])
                    weights.extend(restraints[2])
            return [restraints_vector, jacobian, weights]
        return None


class ScalingRestraints(object):
    """A class to calculate restraints for scaling for an individual dataset, by
    using a single active_parameter_manager."""

    def calculate_jacobian_restraints(self, apm):
        """Calculate jacobian restraints for a single dataset from the scaling model
        components, using a single active_parameter_manager. Return None if no
        restraints, else return a residuals vector and a matrix of size
        n_restrained_parameters x n_parameters, and a weights vector."""
        residual_restraints = self.calculate_restraints(apm)
        if residual_restraints:
            n_restraints = residual_restraints[0].size()
            weights = flex.double([])
            restraints_vector = flex.double([])
            jacobian = sparse.matrix(n_restraints, apm.n_active_params)
            cumul_restr_pos = 0
            for comp in apm.components.values():
                restraints = comp["object"].calculate_jacobian_restraints()
                if restraints:
                    jacobian.assign_block(
                        restraints[1], cumul_restr_pos, comp["start_idx"]
                    )
                    cumul_restr_pos += comp["n_params"]
                    restraints_vector.extend(restraints[0])
                    weights.extend(restraints[2])
            # Return the restraints vector, jacobian and weights (unity as weights
            # contained in individual jacobian/weoghts calculations).
            return [restraints_vector, jacobian, weights]
        return None

    @staticmethod
    def calculate_restraints(apm):
        """Calculate restraints for a single dataset from the scaling model
        components. Return None if no restraints, else return a residuals vector
        length n_restrained_parameters and a gradient vector of length n_parameters
        (of the scaling model for the individual dataset)."""
        residuals = flex.double([])
        gradient_vector = flex.double([])
        for comp in apm.components.values():
            resid = comp["object"].calculate_restraints()
            if resid:
                gradient_vector.extend(resid[1])
                residuals.extend(resid[0])
            else:
                gradient_vector.extend(flex.double(comp["n_params"], 0.0))
        if residuals:
            return [residuals, gradient_vector]
        return None
