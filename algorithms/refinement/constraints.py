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

from __future__ import absolute_import, division
from libtbx.phil import parse
from libtbx.utils import Sorry
from scitbx.array_family import flex
from scitbx import sparse

# PHIL options for constraints
phil_str = '''
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
    .type = ints(value_min=0,size_min=2)

  parameters = None
    .type = strings
    .help = "Constrain specified parameters of each parameterisation by a list"
            "of parameter names to match. Model name prefixes such as"
            "'Detector1' will be ignored as parameterisations are identified"
            "by experiment id"
}

'''

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
    self._constrained_idx = flex.size_t([i for c in self._constrained_gps for i in c])
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

    unconstrained_part = constrained_x[0:self._n_unconstrained_params]
    constrained_part = constrained_x[self._n_unconstrained_params:]

    # update constrained parameter values
    for v, c in zip(constrained_part, self._constraints):
      c.set_constrained_value(v)

    expanded = flex.double([v for c in self._constraints for v in c.get_expanded_values()])

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
      constr_block, 0, self._n_unconstrained_params)

    return constrained_jacobian

class SparseConstraintManager(ConstraintManager):

  def constrain_jacobian(self, jacobian):
    '''sparse matrix version of constrain_jacobian'''

    # select unconstrained columns only
    unconstr_block = jacobian.select_columns(self._unconstrained_idx)

    # create constrained columns
    constr_block = sparse.matrix(jacobian.n_rows, len(self._constraints))

    mask = flex.bool(jacobian.n_rows, True)
    for i, (gp, c) in enumerate(zip(self._constrained_gps, constr_block.cols())):
      # this copies, so c is no longer the matrix column but a new vector
      for j in gp: c += jacobian.col(j)
      # so assign back into the matrix directly
      constr_block[:,i] = c

    # construct the constrained Jacobian
    constrained_jacobian = sparse.matrix(jacobian.n_rows,
      unconstr_block.n_cols + constr_block.n_cols)
    constrained_jacobian.assign_block(unconstr_block, 0, 0)
    constrained_jacobian.assign_block(constr_block, 0, unconstr_block.n_cols)

    return constrained_jacobian
