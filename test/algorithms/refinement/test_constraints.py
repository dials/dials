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

from __future__ import absolute_import, division, print_function
import os
from copy import deepcopy
from libtbx.test_utils import approx_equal
from libtbx import easy_run
from scitbx import sparse
from dials.array_family import flex
from dxtbx.model.experiment_list import ExperimentListFactory
from dials.algorithms.refinement.constraints import EqualShiftConstraint
from dials.algorithms.refinement.constraints import ConstraintManager
from dials.algorithms.refinement.constraints import SparseConstraintManager
from dials.algorithms.refinement.engine import Journal


def test_contraints_manager_simple_test():

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
    mod_constrained_x = constrained_x + 10.0

    # check the expanded parameters are as expected
    expanded = cm.expand_parameters(mod_constrained_x)
    assert x + 10.0 == expanded

    # make a matrix to exercise jacobian compaction
    j = flex.random_double(20 * 10)
    j.reshape(flex.grid(20, 10))

    # for constrained columns, elements that are non-zero in one column are
    # zero in the other columns. Enforce that in this example
    mask2 = flex.bool([True] * 10 + [False] * 10)
    mask4 = ~mask2
    col2 = j.matrix_copy_column(1)
    col2.set_selected(mask2, 0)
    j.matrix_paste_column_in_place(col2, 1)
    col4 = j.matrix_copy_column(3)
    col4.set_selected(mask4, 0)
    j.matrix_paste_column_in_place(col4, 3)

    mask6 = flex.bool([False] * 7 + [True] * 13)
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

    # convert to a sparse matrix to exercise the sparse Jacobian compaction
    j2 = sparse.matrix(20, 10)
    mask = flex.bool(20, True)
    for i, c in enumerate(j2.cols()):
        c.set_selected(mask, j.matrix_copy_column(i))
    assert (j2.as_dense_matrix() == j).all_eq(True)

    cm2 = SparseConstraintManager([c1, c2], len(x))
    cj2 = cm2.constrain_jacobian(j2)

    # ensure dense and sparse calculations give the same result
    assert (cj2.as_dense_matrix() == cj).all_eq(True)

    # construct derivatives of the objective dL/dp from the Jacobian to test
    # constrain_gradient_vector. Here assume unit weights
    dL_dp = [sum(col.as_dense_vector()) for col in j2.cols()]
    constr_dL_dp = cm.constrain_gradient_vector(dL_dp)

    # check constrained values are equal to sums of relevant elements in the
    # original gradient vector
    assert constr_dL_dp[5] == dL_dp[1] + dL_dp[3]
    assert constr_dL_dp[6] == dL_dp[5] + dL_dp[6] + dL_dp[7]


def test_constrained_refinement(dials_regression, run_in_tmpdir):
    """Test joint refinement where two detectors are constrained to enforce a
    differential distance (along the shared initial normal vector) of 1 mm.
    This test can be constructed on the fly from data already in
    dials_regression"""

    # use the 'centroid' data for this test. The 'regularized' experiments are
    # useful because the detector has fast and slow exactly aligned with X, -Y
    # so the distance is exactly along the normal vector and can be altered
    # directly by changing the Z component of the orgin vector
    data_dir = os.path.join(dials_regression, "refinement_test_data", "centroid")
    experiments_path = os.path.join(data_dir, "experiments_XPARM_REGULARIZED.json")
    pickle_path = os.path.join(data_dir, "spot_1000_xds.pickle")

    # load the experiments and spots
    el = ExperimentListFactory.from_json_file(experiments_path, check_format=False)
    rt = flex.reflection_table.from_pickle(pickle_path)

    # adjust the detector distance by -0.5 mm
    detector = el[0].detector
    panel = detector[0]
    fast = panel.get_fast_axis()
    slow = panel.get_slow_axis()
    origin = panel.get_origin()
    panel.set_frame(fast, slow, origin[0:2] + (origin[2] + 0.5,))

    # duplicate the experiment and adjust distance by +1 mm
    e2 = deepcopy(el[0])
    detector = e2.detector
    panel = detector[0]
    fast = panel.get_fast_axis()
    slow = panel.get_slow_axis()
    origin = panel.get_origin()
    panel.set_frame(fast, slow, origin[0:2] + (origin[2] - 1.0,))

    # append to the experiment list and write out
    el.append(e2)
    el.as_json("foo.expt")

    # duplicate the reflections and increment the experiment id
    rt2 = deepcopy(rt)
    rt2["id"] = rt2["id"] + 1

    # concatenate reflections and write out
    rt.extend(rt2)
    rt.as_pickle("foo.refl")

    # set up refinement, constraining the distance parameter
    cmd = (
        "dials.refine foo.expt foo.refl "
        "history=history.json refinement.parameterisation.detector."
        "constraints.parameter=Dist scan_varying=False"
    )
    easy_run.fully_buffered(command=cmd).raise_if_errors()

    # load refinement history
    history = Journal.from_json_file("history.json")
    ref_exp = ExperimentListFactory.from_json_file("refined.expt", check_format=False)

    # we expect 8 steps of constrained refinement
    assert history.get_nrows() == 8

    # get parameter vector from the final step
    pvec = history["parameter_vector"][-1]

    # the constrained parameters have indices 0 and 6 in this case. Check they
    # are still exactly 1 mm apart
    assert pvec[0] == pvec[6] - 1.0

    # NB because the other detector parameters were not also constrained, the
    # refined lab frame distances may not in fact differ by 1 mm. The constraint
    # acts along the initial detector normal vector during composition of a new
    # detector position. After refinement of tilt/twist type rotations,
    # the final distances along the new normal vectors will change
    det1, det2 = ref_exp.detectors()
    p1 = det1[0]
    p2 = det2[0]
    assert approx_equal(p2.get_distance() - p1.get_distance(), 0.9987655)
