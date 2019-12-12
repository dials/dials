from __future__ import absolute_import, division, print_function

import logging

import scitbx.matrix
from cctbx import crystal, sgtbx
from cctbx.crystal_orientation import crystal_orientation
from cctbx.sgtbx import change_of_basis_op, lattice_symmetry, subgroups
from cctbx.sgtbx.bravais_types import bravais_lattice
from rstbx.dps_core.lepage import iotbx_converter
from scitbx.array_family import flex

from dxtbx.model import Crystal

from dials.algorithms.indexing import DialsIndexError

logger = logging.getLogger(__name__)


def metric_supergroup(group):
    return (
        sgtbx.space_group_info(group=group)
        .type()
        .expand_addl_generators_of_euclidean_normalizer(True, True)
        .build_derived_acentric_group()
    )


def find_matching_symmetry(unit_cell, target_space_group, max_delta=5):
    cs = crystal.symmetry(unit_cell=unit_cell, space_group=sgtbx.space_group())
    target_bravais_t = bravais_lattice(
        group=target_space_group.info().reference_setting().group()
    )
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

    # Get list of sub-spacegroups
    subgrs = subgroups.subgroups(lattice_group.info()).groups_parent_setting()

    # Order sub-groups
    sort_values = flex.double()
    for group in subgrs:
        order_z = group.order_z()
        space_group_number = sgtbx.space_group_type(group, False).number()
        assert 1 <= space_group_number <= 230
        sort_values.append(order_z * 1000 + space_group_number)
    perm = flex.sort_permutation(sort_values, reverse=True)

    for i_subgr in perm:
        acentric_subgroup = subgrs[i_subgr]
        acentric_supergroup = metric_supergroup(acentric_subgroup)
        # Make symmetry object: unit-cell + space-group
        # The unit cell is potentially modified to be exactly compatible
        # with the space group symmetry.
        subsym = crystal.symmetry(
            unit_cell=minimum_symmetry.unit_cell(),
            space_group=acentric_subgroup,
            assert_is_compatible_unit_cell=False,
        )
        # Convert subgroup to reference setting
        cb_op_minimum_ref = subsym.space_group_info().type().cb_op()
        ref_subsym = subsym.change_basis(cb_op_minimum_ref)
        # Ignore unwanted groups
        bravais_t = bravais_lattice(group=ref_subsym.space_group())
        if bravais_t != target_bravais_t:
            continue

        # Choose best setting for monoclinic and orthorhombic systems
        cb_op_best_cell = ref_subsym.change_of_basis_op_to_best_cell(
            best_monoclinic_beta=True
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


class SymmetryHandler(object):
    def __init__(self, unit_cell=None, space_group=None, max_delta=5):

        self._max_delta = max_delta
        self.target_symmetry_primitive = None
        self.target_symmetry_reference_setting = None
        self.cb_op_inp_ref = None

        target_unit_cell = unit_cell
        target_space_group = space_group
        if target_space_group is not None:
            target_space_group = target_space_group.build_derived_patterson_group()

        if target_unit_cell is not None or target_space_group is not None:

            if target_unit_cell is not None and target_space_group is not None:
                self._setup_target_unit_cell_and_space_group(
                    target_unit_cell, target_space_group
                )

            elif target_unit_cell is not None:
                self.target_symmetry_reference_setting = crystal.symmetry(
                    unit_cell=target_unit_cell, space_group=sgtbx.space_group()
                )
                self.cb_op_inp_ref = sgtbx.change_of_basis_op()

            elif target_space_group is not None:
                self.cb_op_inp_ref = (
                    target_space_group.info().change_of_basis_op_to_reference_setting()
                )
                self.target_symmetry_reference_setting = crystal.symmetry(
                    space_group=target_space_group.change_basis(self.cb_op_inp_ref)
                )

            self.cb_op_reference_to_primitive = (
                self.target_symmetry_reference_setting.change_of_basis_op_to_primitive_setting()
            )
            if target_unit_cell is not None:
                self.target_symmetry_primitive = self.target_symmetry_reference_setting.change_basis(
                    self.cb_op_reference_to_primitive
                )
            else:
                self.target_symmetry_primitive = crystal.symmetry(
                    space_group=self.target_symmetry_reference_setting.space_group().change_basis(
                        self.cb_op_reference_to_primitive
                    )
                )
            self.cb_op_ref_inp = self.cb_op_inp_ref.inverse()
            self.cb_op_primitive_inp = (
                self.cb_op_ref_inp * self.cb_op_reference_to_primitive.inverse()
            )

            if self.target_symmetry_reference_setting:
                logger.debug(
                    "Target symmetry (reference setting):\n%s",
                    self.target_symmetry_reference_setting,
                )
            if self.target_symmetry_primitive:
                logger.debug(
                    "Target symmetry (primitive cell):\n%s",
                    self.target_symmetry_primitive,
                )
            logger.debug(
                "cb_op reference->primitive: %s", self.cb_op_reference_to_primitive
            )
            logger.debug("cb_op primitive->input: %s", self.cb_op_primitive_inp)

    def _setup_target_unit_cell_and_space_group(
        self, target_unit_cell, target_space_group
    ):

        target_bravais_t = bravais_lattice(
            group=target_space_group.info().reference_setting().group()
        )
        best_subgroup = None
        best_angular_difference = 1e8

        space_groups = [target_space_group]
        if target_space_group.conventional_centring_type_symbol() != "P":
            space_groups.append(sgtbx.space_group())
        for target in space_groups:
            cs = crystal.symmetry(
                unit_cell=target_unit_cell,
                space_group=target,
                assert_is_compatible_unit_cell=False,
            )
            target_best_cell = cs.best_cell().unit_cell()
            subgroups = lattice_symmetry.metric_subgroups(cs, max_delta=0.1)
            for subgroup in subgroups.result_groups:
                bravais_t = bravais_lattice(group=subgroup["ref_subsym"].space_group())
                if bravais_t == target_bravais_t:
                    # allow for the cell to be given as best cell, reference setting
                    # primitive settings, or minimum cell
                    best_subsym = subgroup["best_subsym"]
                    ref_subsym = best_subsym.as_reference_setting()
                    if not (
                        best_subsym.unit_cell().is_similar_to(target_unit_cell)
                        or ref_subsym.unit_cell().is_similar_to(target_unit_cell)
                        or ref_subsym.primitive_setting()
                        .unit_cell()
                        .is_similar_to(target_unit_cell)
                        or best_subsym.primitive_setting()
                        .unit_cell()
                        .is_similar_to(target_unit_cell)
                        or best_subsym.minimum_cell()
                        .unit_cell()
                        .is_similar_to(target_unit_cell.minimum_cell())
                        or best_subsym.unit_cell().is_similar_to(target_best_cell)
                    ):
                        continue
                    if subgroup["max_angular_difference"] < best_angular_difference:
                        best_subgroup = subgroup
                        best_angular_difference = subgroup["max_angular_difference"]

        if best_subgroup is None:
            raise DialsIndexError("Unit cell incompatible with space group")

        cb_op_inp_best = best_subgroup["cb_op_inp_best"]
        best_subsym = best_subgroup["best_subsym"]
        cb_op_best_ref = best_subsym.change_of_basis_op_to_reference_setting()
        self.cb_op_inp_ref = cb_op_best_ref * cb_op_inp_best
        self.target_symmetry_reference_setting = crystal.symmetry(
            unit_cell=target_unit_cell.change_basis(self.cb_op_inp_ref),
            space_group=target_space_group.info().as_reference_setting().group(),
        )

    def apply_symmetry(self, crystal_model):
        if not (
            self.target_symmetry_primitive
            and self.target_symmetry_primitive.space_group()
        ):
            return crystal, sgtbx.change_of_basis_op()

        target_space_group = self.target_symmetry_primitive.space_group()

        A = crystal_model.get_A()

        max_delta = self._max_delta
        items = iotbx_converter(crystal_model.get_unit_cell(), max_delta=max_delta)
        target_sg_ref = target_space_group.info().reference_setting().group()
        best_angular_difference = 1e8
        best_subgroup = None
        for item in items:
            if bravais_lattice(group=target_sg_ref) != bravais_lattice(
                group=item["ref_subsym"].space_group()
            ):
                continue
            if item["max_angular_difference"] < best_angular_difference:
                best_angular_difference = item["max_angular_difference"]
                best_subgroup = item

        if best_subgroup is None:
            return None, None

        cb_op_inp_best = best_subgroup["cb_op_inp_best"]
        orient = crystal_orientation(A, True)
        orient_best = orient.change_basis(
            scitbx.matrix.sqr(cb_op_inp_best.c().as_double_array()[0:9]).transpose()
        )
        constrain_orient = orient_best.constrain(best_subgroup["system"])

        best_subsym = best_subgroup["best_subsym"]
        cb_op_best_ref = best_subsym.change_of_basis_op_to_reference_setting()
        target_sg_best = target_sg_ref.change_basis(cb_op_best_ref.inverse())
        ref_subsym = best_subsym.change_basis(cb_op_best_ref)
        cb_op_ref_primitive = ref_subsym.change_of_basis_op_to_primitive_setting()
        cb_op_best_primitive = cb_op_ref_primitive * cb_op_best_ref
        cb_op_inp_primitive = cb_op_ref_primitive * cb_op_best_ref * cb_op_inp_best

        direct_matrix = constrain_orient.direct_matrix()

        a = scitbx.matrix.col(direct_matrix[:3])
        b = scitbx.matrix.col(direct_matrix[3:6])
        c = scitbx.matrix.col(direct_matrix[6:9])
        model = Crystal(a, b, c, space_group=target_sg_best)
        assert target_sg_best.is_compatible_unit_cell(model.get_unit_cell())

        model = model.change_basis(cb_op_best_primitive)
        return model, cb_op_inp_primitive
