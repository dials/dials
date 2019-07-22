from __future__ import absolute_import, division
from __future__ import print_function

import copy

from cctbx import sgtbx
import scitbx.matrix
from scitbx.math import euler_angles_as_matrix
from dxtbx.model import Crystal

from dials.algorithms.indexing.basis_vector_search import strategies
from dials.algorithms.indexing.basis_vector_search import combinations


def test_combinations(setup_rlp):
    max_cell = 1.3 * max(setup_rlp["crystal_symmetry"].unit_cell().parameters()[:3])
    strategy = strategies.FFT1D(max_cell)
    basis_vectors, used = strategy.find_basis_vectors(setup_rlp["rlp"])

    for target_symmetry in (
        setup_rlp["crystal_symmetry"],
        setup_rlp["crystal_symmetry"]
        .primitive_setting()
        .customized_copy(space_group_info=sgtbx.space_group().info()),
        setup_rlp["crystal_symmetry"]
        .primitive_setting()
        .customized_copy(unit_cell=None),
        setup_rlp["crystal_symmetry"].as_reference_setting(),
    ):

        crystal_models = combinations.candidate_orientation_matrices(
            basis_vectors, max_combinations=50
        )
        filtered_crystal_models = combinations.filter_known_symmetry(
            crystal_models, target_symmetry=target_symmetry
        )
        filtered_crystal_models = list(filtered_crystal_models)

        assert len(filtered_crystal_models)
        for model in filtered_crystal_models:
            if target_symmetry.unit_cell() is not None:
                assert model.get_unit_cell().minimum_cell().is_similar_to(
                    target_symmetry.minimum_cell().unit_cell(),
                    relative_length_tolerance=0.1,
                    absolute_angle_tolerance=5,
                ) or model.get_unit_cell().is_similar_to(
                    target_symmetry.unit_cell(),
                    relative_length_tolerance=0.1,
                    absolute_angle_tolerance=5,
                )
            else:
                target_sg = (
                    target_symmetry.space_group_info().reference_setting().group()
                )
                from cctbx.sgtbx.lattice_symmetry import metric_subgroups

                subgroups = metric_subgroups(
                    model.get_crystal_symmetry(), max_delta=5, bravais_types_only=False
                )
                if not target_sg.build_derived_patterson_group() in [
                    g["ref_subsym"].space_group() for g in subgroups.result_groups
                ]:
                    assert 0


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
