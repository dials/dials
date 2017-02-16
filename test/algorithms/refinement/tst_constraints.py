#!/usr/bin/env cctbx.python

#
#  Copyright (C) (2017) STFC Rutherford Appleton Laboratory, UK.
#
#  Author: David Waterman.
#
#  This code is distributed under the BSD license, a copy of which is
#  included in the root directory of this package.
#

"""
Tests for the constraints system used in refinement

"""

from __future__ import absolute_import, division
from dials.array_family import flex
from dials.algorithms.refinement.constraints import EqualShiftConstraint
from dials.algorithms.refinement.constraints import ConstraintManager
from dials.algorithms.refinement.constraints import SparseConstraintManager
from scitbx import sparse

def test1():

  x = flex.random_double(10)

  # constrain parameters 2 and 4 and 6, 7 and 8
  c1 = EqualShiftConstraint([1, 3], x)
  c2 = EqualShiftConstraint([5, 6, 7], x)

  cm = ConstraintManager([c1, c2], len(x))
  constrained_x = cm.constrain_parameters(x)

  # check the constrained parameters are as expected
  assert len(constrained_x) == 7
  assert constrained_x[5] == flex.mean(x.select([1, 3]))
  assert constrained_x[6] == flex.mean(x[5:8])

  # minimiser would modify the constrained parameters
  mod_constrained_x = constrained_x + 10.

  # check the expanded parameters are as expected
  expanded = cm.expand_parameters(mod_constrained_x)
  assert x + 10. == expanded

  # make a matrix to exercise jacobian compaction
  j = flex.random_double(20*10)
  j.reshape(flex.grid(20,10))

  # for constrained columns, elements that are non-zero in one column are
  # zero in the other columns. Enforce that in this example
  mask2 = flex.bool([True]*10 + [False]*10)
  mask4 = ~mask2
  col2 = j.matrix_copy_column(1)
  col2.set_selected(mask2, 0)
  j.matrix_paste_column_in_place(col2, 1)
  col4 = j.matrix_copy_column(3)
  col4.set_selected(mask4, 0)
  j.matrix_paste_column_in_place(col4, 3)

  mask6 = flex.bool([False]*7 + [True]*13)
  mask7 = mask6.reversed()
  mask8 = ~(mask6 & mask7)
  col6 = j.matrix_copy_column(5)
  col6.set_selected(mask6, 0)
  j.matrix_paste_column_in_place(col6, 5)
  col7 = j.matrix_copy_column(6)
  col7.set_selected(mask7, 0)
  j.matrix_paste_column_in_place(col7, 6)
  col8 = j.matrix_copy_column(7)
  col8.set_selected(mask8, 0)
  j.matrix_paste_column_in_place(col8, 7)

  cj = cm.constrain_jacobian(j)

  # check expected dimensions
  assert cj.all() == (20, 7)

  # check that the constrained columns are equal to sums of the relevant
  # columns in the original Jacobian
  tmp = j.matrix_copy_column(1) + j.matrix_copy_column(3)
  assert (cj.matrix_copy_column(5) == tmp).all_eq(True)

  tmp = j.matrix_copy_column(5) + j.matrix_copy_column(6) + j.matrix_copy_column(7)
  assert (cj.matrix_copy_column(6) == tmp).all_eq(True)

  # convert to a sparse matrix to exercise the sparse jacobian compaction
  j2 = sparse.matrix(20,10)
  mask = flex.bool(20, True)
  for i, c in enumerate(j2.cols()):
    c.set_selected(mask, j.matrix_copy_column(i))
  assert (j2.as_dense_matrix() == j).all_eq(True)

  cm2 = SparseConstraintManager([c1, c2], len(x))
  cj2 = cm2.constrain_jacobian(j2)

  # ensure dense and sparse calculations give the same result
  assert (cj2.as_dense_matrix() == cj).all_eq(True)

  print "OK"

if __name__ == '__main__':

  # simple test of constraint manager
  test1()
