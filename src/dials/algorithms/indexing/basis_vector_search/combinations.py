from __future__ import annotations

import logging
import math
from itertools import combinations as iter_combinations

import numpy as np
from cctbx.sgtbx.bravais_types import bravais_lattice
from cctbx.uctbx.reduction_base import iteration_limit_exceeded
from dxtbx.model import Crystal
from scitbx import matrix as scitbx_matrix

from dials.algorithms.indexing import DialsIndexError
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

    # Build sorted (i<j<k) index array. itertools.combinations already enforces
    # i<j<k, replacing the flex filter. Sort by squared index-norm to match the
    # original ordering (smallest-indexed vectors — best candidates — first).
    idxs = np.array(list(iter_combinations(range(n), 3)), dtype=np.int32)
    norms_sq = (idxs**2).sum(axis=1)
    idxs = idxs[np.argsort(norms_sq, kind="stable")]

    if max_combinations is not None and max_combinations < len(idxs):
        idxs = idxs[:max_combinations]

    half_pi = 0.5 * math.pi
    min_angle = 20 / 180 * math.pi  # 20 degrees, arbitrary cutoff

    # Convert basis vectors to a (n, 3) numpy array.
    bv = np.array([v.elems for v in basis_vectors])

    # Extract all (a, b, c) vector triplets.
    a_arr = bv[idxs[:, 0]]  # (N, 3)
    b_arr = bv[idxs[:, 1]]  # (N, 3)
    c_arr = bv[idxs[:, 2]]  # (N, 3)

    # Precompute squared norms for all vectors in each triplet.
    a_ns = np.einsum("ij,ij->i", a_arr, a_arr)
    b_ns = np.einsum("ij,ij->i", b_arr, b_arr)
    c_ns = np.einsum("ij,ij->i", c_arr, c_arr)

    # Filter 1: angle(a, b) not too close to 0° or 180°.
    ab_dot = np.einsum("ij,ij->i", a_arr, b_arr)
    angle_ab = np.arccos(np.clip(ab_dot / np.sqrt(a_ns * b_ns), -1.0, 1.0))
    m1 = (angle_ab >= min_angle) & ((np.pi - angle_ab) >= min_angle)
    a_arr, b_arr, c_arr = a_arr[m1], b_arr[m1], c_arr[m1]
    b_ns, c_ns, angle_ab = b_ns[m1], c_ns[m1], angle_ab[m1]

    # Flip b (and implicitly a×b) where angle_ab < half_pi so all angles obtuse.
    flip_b = (angle_ab < half_pi)[:, None]
    b_arr = np.where(flip_b, -b_arr, b_arr)

    # Compute a × b after the b flip.
    a_cross_b = np.cross(a_arr, b_arr)  # (M1, 3)
    acb_ns = np.einsum("ij,ij->i", a_cross_b, a_cross_b)

    # Filter 2: angle(a×b, c) not too close to 90° (would give degenerate cell).
    acb_c_dot = np.einsum("ij,ij->i", a_cross_b, c_arr)
    angle_acb_c = np.arccos(np.clip(acb_c_dot / np.sqrt(acb_ns * c_ns), -1.0, 1.0))
    m2 = np.abs(half_pi - angle_acb_c) >= min_angle
    a_arr, b_arr, c_arr = a_arr[m2], b_arr[m2], c_arr[m2]
    a_cross_b, b_ns, c_ns = a_cross_b[m2], b_ns[m2], c_ns[m2]

    # Flip c where angle(b, c) < half_pi so all angles obtuse.
    # alpha is in radians — consistent with half_pi (fixes a deg=True bug in the
    # original code where alpha was computed in degrees but compared to half_pi).
    bc_dot = np.einsum("ij,ij->i", b_arr, c_arr)
    alpha = np.arccos(np.clip(bc_dot / np.sqrt(b_ns * c_ns), -1.0, 1.0))
    flip_c = (alpha < half_pi)[:, None]
    c_arr = np.where(flip_c, -c_arr, c_arr)

    # Ensure right-handed basis: invert all vectors if a×b · c < 0.
    acb_dot_c = np.einsum("ij,ij->i", a_cross_b, c_arr)
    flip_all = (acb_dot_c < 0)[:, None]
    a_arr = np.where(flip_all, -a_arr, a_arr)
    b_arr = np.where(flip_all, -b_arr, b_arr)
    c_arr = np.where(flip_all, -c_arr, c_arr)

    # Crystal creation and Niggli reduction are C++ and cannot be batched.
    # The loop runs only over the fraction of combinations that passed the filters.
    for a_row, b_row, c_row in zip(a_arr, b_arr, c_arr):
        a = scitbx_matrix.col(a_row.tolist())
        b = scitbx_matrix.col(b_row.tolist())
        c = scitbx_matrix.col(c_row.tolist())
        model = Crystal(a, b, c, space_group_symbol="P 1")
        uc = model.get_unit_cell()
        try:
            cb_op_to_niggli = uc.change_of_basis_op_to_niggli_cell()
        except iteration_limit_exceeded as e:
            raise DialsIndexError(e)
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
