from __future__ import annotations

import copy
import logging

from cctbx import sgtbx
from dxtbx.model import ExperimentList
from rstbx.symmetry.constraints import parameter_reduction

from dials.algorithms.scaling.scaling_library import determine_best_unit_cell
from dials.algorithms.symmetry.reindex_to_reference import (
    determine_reindex_operator_against_reference,
)
from dials.array_family import flex
from dials.util.filter_reflections import filtered_arrays_from_experiments_reflections

logger = logging.getLogger(__name__)


def derive_change_of_basis_op(from_hkl, to_hkl):
    # exclude those reflections that we couldn't index
    sel = (to_hkl != (0, 0, 0)) & (from_hkl != (0, 0, 0))
    assert sel.count(True) >= 3  # need minimum of 3 equations ?
    to_hkl = to_hkl.select(sel)
    from_hkl = from_hkl.select(sel)

    # for each miller index, solve a system of linear equations to find the
    # change of basis operator
    h, k, l = to_hkl.as_vec3_double().parts()

    r = []
    from scitbx.lstbx import normal_eqns

    for i in range(3):
        eqns = normal_eqns.linear_ls(3)
        for index, hkl in zip((h, k, l)[i], from_hkl):
            eqns.add_equation(
                right_hand_side=index, design_matrix_row=flex.double(hkl), weight=1
            )
        eqns.solve()
        r.extend(eqns.solution())

    from scitbx import matrix
    from scitbx.math import continued_fraction

    denom = 12
    r = [
        int(denom * continued_fraction.from_real(r_, eps=1e-2).as_rational())
        for r_ in r
    ]
    r = matrix.sqr(r).transpose()
    # print (1/denom)*r

    # now convert into a cctbx change_of_basis_op object
    change_of_basis_op = sgtbx.change_of_basis_op(
        sgtbx.rt_mx(sgtbx.rot_mx(r, denominator=denom))
    ).inverse()
    logger.info(f"discovered change_of_basis_op={change_of_basis_op}")

    # sanity check that this is the right cb_op
    assert (change_of_basis_op.apply(from_hkl) == to_hkl).count(False) == 0

    return change_of_basis_op


def change_of_basis_op_against_reference(
    experiments, reflections, reference_miller_set
):
    # Set some flags to allow filtering, if wanting to reindex against
    # reference with data that has not yet been through integration
    for table in reflections:
        if table.get_flags(table.flags.integrated_sum).count(True) == 0:
            assert "intensity.sum.value" in table, (
                "No 'intensity.sum.value' in reflections"
            )
            table.set_flags(flex.bool(table.size(), True), table.flags.integrated_sum)
    # Make miller array of the datasets
    best_cell = determine_best_unit_cell(experiments)
    filter_logger = logging.getLogger("dials.util.filter_reflections")
    filter_logger.disabled = True
    test_miller_sets = filtered_arrays_from_experiments_reflections(
        experiments, reflections, partiality_threshold=0.4
    )
    test_miller_set = test_miller_sets[0]
    for d in test_miller_sets[1:]:
        test_miller_set = test_miller_set.concatenate(d)
    test_miller_set = test_miller_set.customized_copy(unit_cell=best_cell)
    change_of_basis_op = determine_reindex_operator_against_reference(
        test_miller_set, reference_miller_set
    )
    return change_of_basis_op


def reindex_experiments(
    experiments: ExperimentList,
    cb_op: sgtbx.change_of_basis_op,
    space_group: sgtbx.space_group | None = None,
) -> ExperimentList:
    """
    Reindexes the given experiment list using the provided change of basis operator, and optionally set the space group.

    Args:
        experiments (ExperimentList): The list of experiments to reindex.
        cb_op (sgtbx.change_of_basis_op): The change of basis operator.
        space_group (sgtbx.space_group | None, optional): The space group to set after reindexing. Defaults to None.

    Returns:
        ExperimentList: The reindexed experiments.
    """

    reindexed_experiments = copy.deepcopy(experiments)
    for crystal in reindexed_experiments.crystals():
        cryst_reindexed = copy.deepcopy(crystal)
        if space_group is not None:
            # See also https://github.com/cctbx/cctbx_project/issues/424
            cryst_reindexed.set_space_group(sgtbx.space_group("P 1"))
            cryst_reindexed = cryst_reindexed.change_basis(cb_op)
            cryst_reindexed.set_space_group(space_group)
            S = parameter_reduction.symmetrize_reduce_enlarge(
                cryst_reindexed.get_space_group()
            )
            S.set_orientation(cryst_reindexed.get_B())
            S.symmetrize()
            # Cache the scan-varying A matrices if applicable as these get lost
            # when we call crystal.set_B()
            A_varying = [
                cryst_reindexed.get_A_at_scan_point(i)
                for i in range(cryst_reindexed.num_scan_points)
            ]
            # Update the symmetrized B matrix
            cryst_reindexed.set_B(S.orientation.reciprocal_matrix())
            # Reapply the scan-varying A matrices
            cryst_reindexed.set_A_at_scan_points(A_varying)
        else:
            cryst_reindexed = cryst_reindexed.change_basis(cb_op)
        crystal.update(cryst_reindexed)

    return reindexed_experiments


def reindex_reflections(
    reflections: list[flex.reflection_table],
    change_of_basis_op: sgtbx.change_of_basis_op,
    hkl_offset: list[int] | None = None,
) -> flex.reflection_table:
    """
    Reindexes reflection tables based on a change of basis operation and an optional HKL offset.

    Args:
        reflections (list[flex.reflection_table]): A list of reflection tables to be reindexed.
        change_of_basis_op (sgtbx.change_of_basis_op): The change of basis operation to apply.
        hkl_offset (list[int] | None, optional): An optional HKL offset to apply. Defaults to None.

    Returns:
        flex.reflection_table: The reindexed reflection table.

    Removes reflections whose change of basis results in non-integral indices.
    """

    reflections = flex.reflection_table.concat(reflections)
    miller_indices = reflections["miller_index"]

    if hkl_offset is not None:
        h, k, l = miller_indices.as_vec3_double().parts()
        h += hkl_offset[0]
        k += hkl_offset[1]
        l += hkl_offset[2]
        miller_indices = flex.miller_index(h.iround(), k.iround(), l.iround())
    non_integral_indices = change_of_basis_op.apply_results_in_non_integral_indices(
        miller_indices
    )
    if non_integral_indices.size() > 0:
        logger.info(
            f"Removing {non_integral_indices.size()}/{miller_indices.size()} reflections (change of basis results in non-integral indices)"
        )
    sel = flex.bool(miller_indices.size(), True)
    sel.set_selected(non_integral_indices, False)
    miller_indices_reindexed = change_of_basis_op.apply(miller_indices.select(sel))
    reflections["miller_index"].set_selected(sel, miller_indices_reindexed)
    reflections["miller_index"].set_selected(~sel, (0, 0, 0))

    return reflections
