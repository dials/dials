#!/usr/bin/env python
#
#  constraints.py
#
#  Copyright (C) 2017 Diamond Light Source and STFC Rutherford Appleton
#                     Laboratory, UK.
#
#  Author: David Waterman
#
#  This code is distributed under the BSD license, a copy of which is
#  included in the root directory of this package.

from __future__ import absolute_import, division, print_function

import logging
from functools import reduce

from libtbx.phil import parse
from dials.algorithms.refinement import DialsRefineConfigError
from scitbx import sparse
from scitbx.array_family import flex

logger = logging.getLogger(__name__)

# PHIL options for constraints
phil_str = """
constraints
  .help = "Parameter equal shift constraints to use in refinement."
  .expert_level = 2
  .multiple = True
{
  id = None
    .help = "Select only the specified experiments when looking up which"
            "parameterisations to apply the constraint to. If an identified"
            "parameterisation affects multiple experiments then the index"
            "of any one of those experiments suffices to identify that"
            "parameterisation. If None (the default) then constraints will be"
            "applied to all parameterisations of this type."
    .type = ints(value_min=0)

  parameter = None
    .type = str
    .help = "Identify which parameter of each parameterisation to constrain by"
            "a (partial) parameter name to match. Model name prefixes such as"
            "'Detector1' will be ignored as parameterisations are already"
            "identified by experiment id"
}

"""

phil_scope = parse(phil_str)


class EqualShiftConstraint(object):
    """A single constraint between parameters of the same type in different
    parameterisations"""

    def __init__(self, indices, parameter_vector):

        self.indices = indices
        parameter_vector = flex.double(parameter_vector)
        self.constrained_value = flex.mean(parameter_vector.select(indices))
        self._shifts = parameter_vector.select(indices) - self.constrained_value

    def set_constrained_value(self, val):
        self.constrained_value = val

    def get_expanded_values(self):
        return self.constrained_value + self._shifts


class ConstraintManager(object):
    def __init__(self, constraints, n_full_params):

        self._constraints = constraints

        # constraints should be a list of EqualShiftConstraint objects
        assert len(self._constraints) > 0

        self._n_full_params = n_full_params
        full_idx = flex.size_t_range(n_full_params)
        self._constrained_gps = [c.indices for c in self._constraints]
        self._constrained_idx = flex.size_t(
            [i for c in self._constrained_gps for i in c]
        )
        keep = flex.bool(self._n_full_params, True)
        keep.set_selected(self._constrained_idx, False)
        self._unconstrained_idx = full_idx.select(keep)
        self._n_unconstrained_params = len(self._unconstrained_idx)

    def constrain_parameters(self, x):

        assert len(x) == self._n_full_params

        constrained_vals = flex.double([c.constrained_value for c in self._constraints])

        # select unconstrained parameters only
        unconstrained_x = x.select(self._unconstrained_idx)

        constrained_x = flex.double.concatenate(unconstrained_x, constrained_vals)

        return constrained_x

    def expand_parameters(self, constrained_x):

        unconstrained_part = constrained_x[0 : self._n_unconstrained_params]
        constrained_part = constrained_x[self._n_unconstrained_params :]

        # update constrained parameter values
        for v, c in zip(constrained_part, self._constraints):
            c.set_constrained_value(v)

        expanded = flex.double(
            [v for c in self._constraints for v in c.get_expanded_values()]
        )

        full_x = flex.double(self._n_full_params)
        full_x.set_selected(self._unconstrained_idx, unconstrained_part)
        full_x.set_selected(self._constrained_idx, expanded)

        return full_x

    def constrain_jacobian(self, jacobian):

        # set up result matrix
        nrow = jacobian.all()[0]
        ncol = self._n_unconstrained_params + len(self._constraints)
        constrained_jacobian = flex.double(flex.grid(nrow, ncol))

        # create constrained columns
        constr_block = flex.double(flex.grid(nrow, len(self._constraints)))
        for i, gp in enumerate(self._constrained_gps):
            cols = [jacobian.matrix_copy_column(j) for j in gp]
            sum_col = reduce(lambda x, y: x + y, cols)
            constr_block.matrix_paste_column_in_place(sum_col, i)

        # copy unconstrained columns into the result
        for i, j in enumerate(self._unconstrained_idx):
            col = jacobian.matrix_copy_column(j)
            constrained_jacobian.matrix_paste_column_in_place(col, i)

        # copy the constrained block into the result
        constrained_jacobian.matrix_paste_block_in_place(
            constr_block, 0, self._n_unconstrained_params
        )

        return constrained_jacobian

    def constrain_gradient_vector(self, grad):

        # extract unconstrained gradients into the result
        result = []
        result = list(flex.double(grad).select(self._unconstrained_idx))

        # append constrained gradients
        for i, gp in enumerate(self._constrained_gps):
            vals = [grad[j] for j in gp]
            result.append(sum(vals))

        return result


class SparseConstraintManager(ConstraintManager):
    def constrain_jacobian(self, jacobian):
        """sparse matrix version of constrain_jacobian"""

        # select unconstrained columns only
        unconstr_block = jacobian.select_columns(self._unconstrained_idx)

        # create constrained columns
        constr_block = sparse.matrix(jacobian.n_rows, len(self._constraints))

        for i, (gp, c) in enumerate(zip(self._constrained_gps, constr_block.cols())):
            # this copies, so c is no longer the matrix column but a new vector
            for j in gp:
                c += jacobian.col(j)
            # so assign back into the matrix directly
            constr_block[:, i] = c

        # construct the constrained Jacobian
        constrained_jacobian = sparse.matrix(
            jacobian.n_rows, unconstr_block.n_cols + constr_block.n_cols
        )
        constrained_jacobian.assign_block(unconstr_block, 0, 0)
        constrained_jacobian.assign_block(constr_block, 0, unconstr_block.n_cols)

        return constrained_jacobian


class ConstraintManagerFactory(object):
    """Build equal shift constraints as requested in params and package into
    a constraints manager to be linked to the Refinery"""

    def __init__(self, refinement_phil, pred_param, sparse=False):

        self._params = refinement_phil
        self._pred_param = pred_param

        # full parameter names and values
        self._all_names = self._pred_param.get_param_names()
        self._all_vals = self._pred_param.get_param_vals()

        return

    def build_constraint(self, constraint_scope, parameterisation, model_type):
        """Create a constraint for a single parameter specified by
        constraint_scope"""

        if constraint_scope.id is None:
            # get one experiment id for each parameterisation to apply to all
            constraint_scope.id = [e.get_experiment_ids()[0] for e in parameterisation]

        # find which parameterisations are involved, and if any are scan-varying
        # how many sample points there are
        prefixes = []
        n_samples = 0
        for i, p in enumerate(parameterisation):
            if hasattr(p, "num_samples"):
                ns = p.num_samples()
                if n_samples == 0:
                    n_samples = ns
                if ns != n_samples:
                    raise DialsRefineConfigError(
                        "Constraints cannot be created between scan-varying "
                        "parameterisations when these have a different number of "
                        "sample points."
                    )
            for j in p.get_experiment_ids():
                if j in constraint_scope.id:
                    prefixes.append(model_type + "{}".format(j + 1))
                    break

        # ignore model name prefixes
        import re

        patt1 = re.compile("^" + model_type + "[0-9]+")
        pname = patt1.sub("", constraint_scope.parameter)

        # Use a regex to find the parameters to constrain from a list of all the
        # parameter names. There are multiple parts to this. The first part
        # identifies the relevant model type and parameterisation ordinal index,
        # accepting those that were chosen according to the supplied experiment
        # ids. The next part allows for additional text, like 'Group1' that may
        # be used by a multi-panel detector parameterisation. Then the parameter
        # name itself, like 'Dist'. Finally, to accommodate scan-varying
        # parameterisations, suffixes like '_sample0' and '_sample1' are
        # distinguished so that these are constrained separately.
        for i in range(max(n_samples, 1)):
            patt2 = re.compile(
                "^("
                + "|".join(prefixes)
                + r"){1}(?![0-9])(\w*"
                + pname
                + ")(_sample{})?$".format(i)
            )
            indices = [j for j, s in enumerate(self._all_names) if patt2.match(s)]
            if len(indices) == 1:
                continue
            logger.debug(
                "\nThe following parameters will be constrained "
                "to enforce equal shifts at each step of refinement:"
            )
            for k in indices:
                logger.debug(self._all_names[k])
        return EqualShiftConstraint(indices, self._all_vals)

    def __call__(self):

        # shorten options path
        options = self._params.refinement.parameterisation

        # shorten options paths further for individual parameterisation types
        detector_c = options.detector.constraints
        beam_c = options.beam.constraints
        orientation_c = options.crystal.orientation.constraints
        cell_c = options.crystal.unit_cell.constraints

        # quit early if there are no constraints to apply
        n_constraints = sum(
            [len(e) for e in [detector_c, beam_c, orientation_c, cell_c]]
        )
        if n_constraints == 0:
            return None

        logger.debug("\nConfiguring constraints")

        # list of constraint objects to build up
        constraints = []

        detector_p = self._pred_param.get_detector_parameterisations()
        beam_p = self._pred_param.get_beam_parameterisations()
        orientation_p = self._pred_param.get_crystal_orientation_parameterisations()
        cell_p = self._pred_param.get_crystal_unit_cell_parameterisations()

        for constr in detector_c:
            constraints.append(self.build_constraint(constr, detector_p, "Detector"))
        for constr in beam_c:
            constraints.append(self.build_constraint(constr, beam_p, "Beam"))
        for constr in orientation_c:
            constraints.append(self.build_constraint(constr, orientation_p, "Crystal"))
        for constr in cell_c:
            constraints.append(self.build_constraint(constr, cell_p, "Crystal"))

        if len(constraints) == 0:
            return None

        # return constraints manager
        if options.sparse:
            return SparseConstraintManager(constraints, len(self._all_vals))
        else:
            return ConstraintManager(constraints, len(self._all_vals))
