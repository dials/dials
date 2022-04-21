from __future__ import annotations

import copy
import logging

import scitbx.matrix
from cctbx import crystal, sgtbx, uctbx
from cctbx.sgtbx.lattice_symmetry import metric_subgroups
from dxtbx.model import Crystal
from scitbx.math import euler_angles_as_matrix

from dials.algorithms.indexing.basis_vector_search import FFT1D, combinations
from dials.algorithms.indexing.symmetry import find_matching_symmetry


def test_combinations(setup_rlp):
    max_cell = 1.3 * max(setup_rlp["crystal_symmetry"].unit_cell().parameters()[:3])
    strategy = FFT1D(max_cell)
    basis_vectors, used = strategy.find_basis_vectors(setup_rlp["rlp"])

    target_symmetry_primitive = (
        setup_rlp["crystal_symmetry"]
        .primitive_setting()
        .customized_copy(space_group_info=sgtbx.space_group().info())
    )
    target_symmetry_sg_only = (
        setup_rlp["crystal_symmetry"]
        .primitive_setting()
        .customized_copy(unit_cell=None)
    )
    target_symmetry_ref = setup_rlp["crystal_symmetry"].as_reference_setting()

    for target_symmetry in (
        setup_rlp["crystal_symmetry"],
        target_symmetry_primitive,
        target_symmetry_sg_only,
        target_symmetry_ref,
    ):

        crystal_models = combinations.candidate_orientation_matrices(
            basis_vectors, max_combinations=50
        )
        crystal_models = list(crystal_models)
        filtered_crystal_models = combinations.filter_known_symmetry(
            crystal_models, target_symmetry=target_symmetry
        )
        filtered_crystal_models = list(filtered_crystal_models)

        assert filtered_crystal_models
        for model in filtered_crystal_models:
            best_subgroup = find_matching_symmetry(
                model.get_unit_cell(), target_symmetry.space_group()
            )
            if target_symmetry.unit_cell() is not None:
                assert best_subgroup["best_subsym"].unit_cell().is_similar_to(
                    setup_rlp["crystal_symmetry"]
                    .as_reference_setting()
                    .best_cell()
                    .unit_cell(),
                    relative_length_tolerance=0.1,
                    absolute_angle_tolerance=5,
                ) or best_subgroup[
                    "best_subsym"
                ].minimum_cell().unit_cell().is_similar_to(
                    setup_rlp["crystal_symmetry"]
                    .as_reference_setting()
                    .best_cell()
                    .minimum_cell()
                    .unit_cell(),
                    relative_length_tolerance=0.1,
                    absolute_angle_tolerance=5,
                )
            else:
                target_sg = (
                    target_symmetry.space_group_info().reference_setting().group()
                )

                subgroups = metric_subgroups(
                    model.get_crystal_symmetry(), max_delta=5, bravais_types_only=False
                )
                assert target_sg.build_derived_patterson_group() in [
                    g["ref_subsym"].space_group() for g in subgroups.result_groups
                ]


def test_filter_known_symmetry_no_matches(caplog):
    caplog.set_level(logging.DEBUG)
    unit_cell = uctbx.unit_cell((10, 10, 10, 90, 90, 90))
    crystal_models = []
    # the reciprocal matrix
    B = scitbx.matrix.sqr(unit_cell.fractionalization_matrix()).transpose()
    for _ in range(10):
        crystal_models.append(
            Crystal(B, space_group=sgtbx.space_group(), reciprocal=True)
        )

    target_symmetry = crystal.symmetry(
        unit_cell=(15, 15, 15, 85, 85, 85), space_group=sgtbx.space_group()
    )
    list(
        combinations.filter_known_symmetry(
            crystal_models, target_symmetry=target_symmetry
        )
    )
    assert "Rejecting crystal model inconsistent with input symmetry" in caplog.text
    assert (
        "No crystal models remaining after comparing with known symmetry" in caplog.text
    )


def test_filter_similar_orientations():
    space_group = sgtbx.space_group()
    unit_cell = space_group.info().any_compatible_unit_cell(volume=1000)

    crystal_models = []
    # the reciprocal matrix
    B = scitbx.matrix.sqr(unit_cell.fractionalization_matrix()).transpose()
    for i in range(10):
        U = euler_angles_as_matrix((i * 10, 0, 0), deg=True)
        crystal_models.append(Crystal(U * B, space_group=space_group, reciprocal=True))

    other_crystal_models = crystal_models[5:]
    crystal_models = crystal_models[:5]

    # there should be no similar orientations to be filtered out
    filtered = list(
        combinations.filter_similar_orientations(crystal_models, other_crystal_models)
    )
    assert filtered == crystal_models

    # now add similar crystal model rotated by 3-degrees
    cryst = copy.deepcopy(other_crystal_models[0])
    cryst.set_U(
        scitbx.matrix.sqr(cryst.get_U()) * euler_angles_as_matrix((3, 0, 0), deg=True)
    )
    crystal_models.append(cryst)
    filtered = combinations.filter_similar_orientations(
        crystal_models, other_crystal_models
    )
    assert list(filtered) == crystal_models[:5]

    # with lower minimum_angular_separation no models should be filtered out
    filtered = combinations.filter_similar_orientations(
        crystal_models, other_crystal_models, minimum_angular_separation=2
    )
    assert list(filtered) == crystal_models
