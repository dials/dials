from __future__ import annotations

import logging
import math

from cctbx.sgtbx.bravais_types import bravais_lattice
from dxtbx.model import Crystal
from scitbx.array_family import flex

from dials.algorithms.indexing.compare_orientation_matrices import (
    difference_rotation_matrix_axis_angle,
)
from dials.algorithms.indexing.symmetry import find_matching_symmetry

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
            unit cell comparison (default value is 0.1).
        absolute_angle_tolerance (float): Angular tolerance (in degrees) in unit cell
            comparison (default value is 5).
        max_delta (float): Maximum allowed Le Page delta used in searching for basis
            vector combinations that are consistent with the given symmetry (default
            value is 5).
    """

    n_matched = 0

    cb_op_ref_to_primitive = target_symmetry.change_of_basis_op_to_primitive_setting()

    if target_symmetry.unit_cell() is not None:
        target_symmetry_primitive = target_symmetry.change_basis(cb_op_ref_to_primitive)
    else:
        target_symmetry_primitive = target_symmetry.customized_copy(
            space_group_info=target_symmetry.space_group_info().change_basis(
                cb_op_ref_to_primitive
            )
        )
    target_bravais_str = str(
        bravais_lattice(
            group=target_symmetry_primitive.space_group_info()
            .reference_setting()
            .group()
        )
    )

    for model in crystal_models:
        uc = model.get_unit_cell()
        best_subgroup = find_matching_symmetry(
            uc, None, max_delta=max_delta, target_bravais_str=target_bravais_str
        )
        if best_subgroup is not None:
            if target_symmetry.unit_cell() is not None and not (
                best_subgroup["best_subsym"]
                .unit_cell()
                .is_similar_to(
                    target_symmetry.as_reference_setting().best_cell().unit_cell(),
                    relative_length_tolerance=relative_length_tolerance,
                    absolute_angle_tolerance=absolute_angle_tolerance,
                )
            ):
                logger.debug(
                    "Rejecting crystal model inconsistent with input symmetry:\n"
                    f"  Unit cell: {str(model.get_unit_cell())}"
                )
                continue

            n_matched += 1
            yield model
    if not n_matched:
        logger.warning(
            "No crystal models remaining after comparing with known symmetry"
        )


def filter_similar_orientations(
    crystal_models, other_crystal_models, minimum_angular_separation=5
):
    for cryst in crystal_models:
        orientation_too_similar = False
        for cryst_a in other_crystal_models:
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
