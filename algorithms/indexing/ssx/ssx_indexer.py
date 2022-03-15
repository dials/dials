from __future__ import annotations

import logging

from cctbx import crystal, sgtbx
from cctbx.sgtbx import change_of_basis_op, subgroups
from cctbx.sgtbx.bravais_types import bravais_lattice
from scitbx.array_family import flex

from dials.algorithms.indexing.basis_vector_search import combinations
from dials.algorithms.indexing.lattice_search import BasisVectorSearch
from dials.algorithms.indexing.stills_indexer import StillsIndexer
from dials.algorithms.indexing.symmetry import metric_supergroup

logger = logging.getLogger(__name__)


def calc_acentric_subgroups(lattice_group_info, target_bravais_t):
    subgrs = subgroups.subgroups(lattice_group_info).groups_parent_setting()

    sort_values = flex.double()
    for group in subgrs:
        order_z = group.order_z()
        space_group_number = (
            group.info().type().number()
        )  # sgtbx.space_group_type(group, False).number()
        assert 1 <= space_group_number <= 230
        sort_values.append(order_z * 1000 + space_group_number)
    perm = flex.sort_permutation(sort_values, reverse=True)

    good_acentric_subgroups = []
    good_acentric_supergroups = []
    cb_ops = []

    for i_perm in perm:
        ac_sg = subgrs[i_perm]
        cb_op_minimum_ref = ac_sg.info().type().cb_op()

        ac_sg_ref = ac_sg.change_basis(cb_op_minimum_ref)
        bravais_t = bravais_lattice(group=ac_sg_ref)
        if bravais_t != target_bravais_t:
            continue
        good_acentric_subgroups.append(ac_sg)
        good_acentric_supergroups.append(metric_supergroup(ac_sg))
        cb_ops.append(cb_op_minimum_ref)

    return (good_acentric_subgroups, good_acentric_supergroups, cb_ops)


class AcentricSubgroupsForLatticeGroup:

    instance = None

    class _AcentricSubgroupsForLatticeGroup:
        def __init__(self, target_bravais_t):
            self._cached_results = {}
            self._target_bravais_t = target_bravais_t

        def get_data(self, info_str):
            try:
                return self._cached_results[info_str]
            except KeyError:
                return None

        def set_data(self, info_str, data):
            self._cached_results[info_str] = data

    def __new__(cls, target_bravais_t):
        if not AcentricSubgroupsForLatticeGroup.instance:
            AcentricSubgroupsForLatticeGroup.instance = (
                AcentricSubgroupsForLatticeGroup._AcentricSubgroupsForLatticeGroup(
                    target_bravais_t
                )
            )
        return AcentricSubgroupsForLatticeGroup.instance

    def __getattr__(self, name):
        return getattr(self.instance, name)

    def __setattr__(self, name):
        return setattr(self.instance, name)


def find_matching_symmetry(
    unit_cell,
    target_bravais_t,
    max_delta=5,
    best_monoclinic_beta=True,
):
    cs = crystal.symmetry(unit_cell=unit_cell, space_group=sgtbx.space_group())

    best_subgroup = None
    best_angular_difference = 1e8

    # code based on cctbx/sgtbx/lattice_symmetry.py but optimised to only
    # look at subgroups with the correct bravais type

    input_symmetry = cs
    # Get cell reduction operator
    cb_op_inp_minimum = input_symmetry.change_of_basis_op_to_minimum_cell()

    # New symmetry object with changed basis
    minimum_symmetry = input_symmetry.change_basis(cb_op_inp_minimum)

    # Get highest symmetry compatible with lattice
    lattice_group = sgtbx.lattice_symmetry_group(
        minimum_symmetry.unit_cell(),
        max_delta=max_delta,
        enforce_max_delta_for_generated_two_folds=True,
    )
    lattice_group_info = lattice_group.info()
    info_str = lattice_group_info.type().lookup_symbol()

    # Get list of sub-spacegroups
    subgroup_cacher = AcentricSubgroupsForLatticeGroup(target_bravais_t)
    cached_data = subgroup_cacher.get_data(info_str)
    if not cached_data:
        a, b, c = calc_acentric_subgroups(lattice_group_info, target_bravais_t)
        subgroup_cacher.set_data(info_str, (a, b, c))
    else:
        a, b, c = cached_data
    good_acentric_subgroups, good_acentric_supergroups, cb_ops = (a, b, c)

    for acentric_subgroup, acentric_supergroup, cb_op in zip(
        good_acentric_subgroups,
        good_acentric_supergroups,
        cb_ops,
    ):
        # Convert subgroup to reference setting
        cb_op_minimum_ref = cb_op
        # Make symmetry object: unit-cell + space-group
        # The unit cell is potentially modified to be exactly compatible
        # with the space group symmetry.
        subsym = crystal.symmetry(
            unit_cell=minimum_symmetry.unit_cell(),
            space_group=acentric_subgroup,
            assert_is_compatible_unit_cell=False,
        )
        ref_subsym = subsym.change_basis(cb_op_minimum_ref)

        # Choose best setting for monoclinic and orthorhombic systems
        cb_op_best_cell = ref_subsym.change_of_basis_op_to_best_cell(
            best_monoclinic_beta=best_monoclinic_beta
        )

        best_subsym = ref_subsym.change_basis(cb_op_best_cell)
        # Total basis transformation
        cb_op_best_cell = change_of_basis_op(
            str(cb_op_best_cell), stop_chars="", r_den=144, t_den=144
        )
        cb_op_minimum_ref = change_of_basis_op(
            str(cb_op_minimum_ref), stop_chars="", r_den=144, t_den=144
        )
        cb_op_inp_minimum = change_of_basis_op(
            str(cb_op_inp_minimum), stop_chars="", r_den=144, t_den=144
        )
        cb_op_inp_best = cb_op_best_cell * cb_op_minimum_ref * cb_op_inp_minimum
        # Use identity change-of-basis operator if possible
        if best_subsym.unit_cell().is_similar_to(input_symmetry.unit_cell()):
            cb_op_corr = cb_op_inp_best.inverse()
            try:
                best_subsym_corr = best_subsym.change_basis(cb_op_corr)
            except RuntimeError as e:
                if str(e).find("Unsuitable value for rational rotation matrix.") < 0:
                    raise
            else:
                if best_subsym_corr.space_group() == best_subsym.space_group():
                    cb_op_inp_best = cb_op_corr * cb_op_inp_best

        max_angular_difference = sgtbx.lattice_symmetry_find_max_delta(
            reduced_cell=minimum_symmetry.unit_cell(), space_group=acentric_supergroup
        )

        if max_angular_difference < best_angular_difference:
            best_angular_difference = max_angular_difference
            best_subgroup = {
                "subsym": subsym,
                "ref_subsym": ref_subsym,
                "best_subsym": best_subsym,
                "cb_op_inp_best": cb_op_inp_best,
                "max_angular_difference": max_angular_difference,
            }

    if best_subgroup is not None:
        return best_subgroup


def filter_known_symmetry(
    symmetry_handler,
    crystal_models,
    relative_length_tolerance,
    absolute_angle_tolerance,
    max_delta,
):

    n_matched = 0
    for model in crystal_models:
        uc = model.get_unit_cell()
        best_subgroup = find_matching_symmetry(
            uc,
            symmetry_handler.target_bravais_t,
            max_delta=max_delta,
        )
        if best_subgroup is not None:
            if symmetry_handler.unit_cell is not None and not (
                best_subgroup["best_subsym"]
                .unit_cell()
                .is_similar_to(
                    symmetry_handler.reference_best_cell,
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


class SSXIndexerBasisVectorSearch(StillsIndexer, BasisVectorSearch):
    def find_candidate_orientation_matrices(self, candidate_basis_vectors):

        candidate_crystal_models = combinations.candidate_orientation_matrices(
            candidate_basis_vectors,
            max_combinations=self.params.basis_vector_combinations.max_combinations,
        )

        if self._symmetry_handler.target_symmetry_reference_setting is not None:
            target_symmetry = self._symmetry_handler.target_symmetry_reference_setting
        elif self._symmetry_handler.target_symmetry_primitive is not None:
            target_symmetry = self._symmetry_handler.target_symmetry_primitive
        else:
            target_symmetry = None
        if target_symmetry is not None:
            candidate_crystal_models = filter_known_symmetry(
                self._symmetry_handler,
                candidate_crystal_models,
                relative_length_tolerance=self.params.known_symmetry.relative_length_tolerance,
                absolute_angle_tolerance=self.params.known_symmetry.absolute_angle_tolerance,
                max_delta=self.params.known_symmetry.max_delta,
            )

        if self.refined_experiments is not None and len(self.refined_experiments) > 0:
            candidate_crystal_models = combinations.filter_similar_orientations(
                candidate_crystal_models,
                self.refined_experiments.crystals(),
                minimum_angular_separation=self.params.multiple_lattice_search.minimum_angular_separation,
            )

        return candidate_crystal_models
