from __future__ import absolute_import, division
from __future__ import print_function
import logging

import math

from scitbx.array_family import flex
from cctbx import sgtbx
from rstbx.indexing_api import tools
from dxtbx.model import Crystal

from dials.algorithms.indexing.symmetry import find_matching_symmetry
from dials.algorithms.indexing.compare_orientation_matrices import (
    difference_rotation_matrix_axis_angle,
)

logger = logging.getLogger(__name__)


def candidate_orientation_matrices(basis_vectors, max_combinations=None):

    # select unique combinations of input vectors to test
    # the order of combinations is such that combinations comprising vectors
    # nearer the beginning of the input list will appear before combinations
    # comprising vectors towards the end of the list
    n = len(basis_vectors)
    # hardcoded limit on number of vectors, fixes issue #72
    # https://github.com/dials/dials/issues/72
    n = min(n, 100)
    basis_vectors = basis_vectors[:n]
    combinations = flex.vec3_int(flex.nested_loop((n, n, n)))
    combinations = combinations.select(
        flex.sort_permutation(combinations.as_vec3_double().norms())
    )

    # select only those combinations where j > i and k > j
    i, j, k = combinations.as_vec3_double().parts()
    sel = flex.bool(len(combinations), True)
    sel &= j > i
    sel &= k > j
    combinations = combinations.select(sel)

    if max_combinations is not None and max_combinations < len(combinations):
        combinations = combinations[:max_combinations]

    half_pi = 0.5 * math.pi
    min_angle = 20 / 180 * math.pi  # 20 degrees, arbitrary cutoff
    for i, j, k in combinations:
        a = basis_vectors[i]
        b = basis_vectors[j]
        angle = a.angle(b)
        if angle < min_angle or (math.pi - angle) < min_angle:
            continue
        a_cross_b = a.cross(b)
        gamma = a.angle(b)
        if gamma < half_pi:
            # all angles obtuse if possible please
            b = -b
            a_cross_b = -a_cross_b
        c = basis_vectors[k]
        if abs(half_pi - a_cross_b.angle(c)) < min_angle:
            continue
        alpha = b.angle(c, deg=True)
        if alpha < half_pi:
            c = -c
        if a_cross_b.dot(c) < 0:
            # we want right-handed basis set, therefore invert all vectors
            a = -a
            b = -b
            c = -c
        model = Crystal(a, b, c, space_group_symbol="P 1")
        uc = model.get_unit_cell()
        cb_op_to_niggli = uc.change_of_basis_op_to_niggli_cell()
        model = model.change_basis(cb_op_to_niggli)

        uc = model.get_unit_cell()
        params = uc.parameters()
        if uc.volume() > (params[0] * params[1] * params[2] / 100):
            # unit cell volume cutoff from labelit 2004 paper
            yield model


def filter_known_symmetry(
    crystal_models,
    target_symmetry,
    relative_length_tolerance=0.1,
    absolute_angle_tolerance=5,
    max_delta=5,
):
    """Filter crystal models for known symmetry.

    Args:
        crystal_models (list): A list of :class:`dxtbx.model.Crystal` objects.
        target_symmetry (cctbx.crystal.symmetry): The target symmetry for filtering.
        relative_length_tolerance (float): Relative tolerance for unit cell lengths in
            unit cell comparision (default value is 0.1).
        absolute_angle_tolerance (float): Angular tolerance (in degrees) in unit cell
            comparison (default value is 5).
        max_delta (float): Maximum allowed Le Page delta used in searching for basis
            vector combinations that are consistent with the given symmetry (default
            value is 5).
    """

    cb_op_ref_to_primitive = target_symmetry.change_of_basis_op_to_primitive_setting()

    if target_symmetry.unit_cell() is not None:
        target_symmetry_primitive = target_symmetry.change_basis(cb_op_ref_to_primitive)
        target_min_cell = target_symmetry_primitive.minimum_cell().unit_cell()
    else:
        target_symmetry_primitive = target_symmetry.customized_copy(
            space_group_info=target_symmetry.space_group_info().change_basis(
                cb_op_ref_to_primitive
            )
        )
        target_min_cell = None

    transformations = [
        sgtbx.change_of_basis_op(
            str(sgtbx.rt_mx(sgtbx.rot_mx(T["trans"].transpose().as_int())))
        )
        for T in tools.R
        if T["mod"] < 5
    ]
    transformations = []  # XXX temporarily disable cell doubling checks
    transformations.insert(0, sgtbx.change_of_basis_op())

    for model in crystal_models:
        best_model = None
        for T in transformations:
            uc = model.get_unit_cell()
            det = T.c().r().determinant()
            if target_symmetry_primitive.unit_cell() is None:
                if det > 1:
                    break
            else:
                primitive_volume = target_symmetry_primitive.unit_cell().volume()
            if det > 1 and abs(uc.volume() / primitive_volume - det) < 1e-1:
                uc = uc.change_basis(T)
            best_subgroup = find_matching_symmetry(
                uc, target_symmetry_primitive.space_group(), max_delta=max_delta
            )
            cb_op_extra = None
            if best_subgroup is None:
                if not cb_op_ref_to_primitive.is_identity_op():
                    # if we have been told we have a centred unit cell check that
                    # indexing hasn't found the centred unit cell instead of the
                    # primitive cell
                    best_subgroup = find_matching_symmetry(
                        uc,
                        target_symmetry.space_group().build_derived_point_group(),
                        max_delta=max_delta,
                    )
                    cb_op_extra = cb_op_ref_to_primitive
                    if best_subgroup is None:
                        continue
                else:
                    continue
            cb_op_inp_best = best_subgroup["cb_op_inp_best"]
            best_subsym = best_subgroup["best_subsym"]
            cb_op_best_ref = best_subsym.change_of_basis_op_to_reference_setting()
            ref_subsym = best_subsym.change_basis(cb_op_best_ref)
            cb_op_ref_primitive = ref_subsym.change_of_basis_op_to_primitive_setting()
            cb_op_to_primitive = (
                cb_op_ref_primitive * cb_op_best_ref * cb_op_inp_best * T
            )
            if cb_op_extra is not None:
                cb_op_to_primitive = cb_op_extra * cb_op_to_primitive
            best_model = model.change_basis(cb_op_to_primitive)

            if target_symmetry_primitive.unit_cell() is not None and not (
                best_model.get_unit_cell().is_similar_to(
                    target_symmetry_primitive.unit_cell(),
                    relative_length_tolerance=relative_length_tolerance,
                    absolute_angle_tolerance=absolute_angle_tolerance,
                )
                or best_model.get_unit_cell()
                .minimum_cell()
                .is_similar_to(
                    target_min_cell,
                    relative_length_tolerance=relative_length_tolerance,
                    absolute_angle_tolerance=absolute_angle_tolerance,
                )
            ):
                best_model = None
                continue
            else:
                break

        if best_model is not None:
            yield best_model


def filter_similar_orientations(
    crystal_models, other_crystal_models, minimum_angular_separation=5
):
    for cryst in crystal_models:
        orientation_too_similar = False
        for i_a, cryst_a in enumerate(other_crystal_models):
            R_ab, axis, angle, cb_op_ab = difference_rotation_matrix_axis_angle(
                cryst_a, cryst
            )
            if abs(angle) < minimum_angular_separation:  # degrees
                orientation_too_similar = True
                break
        if orientation_too_similar:
            logger.debug("skipping crystal: too similar to other crystals")
            continue
        yield cryst
