from __future__ import annotations

import logging

import scitbx.matrix
from cctbx import crystal, sgtbx
from cctbx.crystal_orientation import crystal_orientation
from cctbx.sgtbx import change_of_basis_op, subgroups
from cctbx.sgtbx.bravais_types import bravais_lattice
from dxtbx.model import Crystal
from rstbx.dps_core.lepage import iotbx_converter
from scitbx.array_family import flex

logger = logging.getLogger(__name__)


def metric_supergroup(group):
    return (
        sgtbx.space_group_info(group=group)
        .type()
        .expand_addl_generators_of_euclidean_normalizer(True, True)
        .build_derived_acentric_group()
    )


def groups_cache(fn):
    class MultiClassCache(object):

        "A set of caches for different bravais types"

        instances = {}

        def __new__(cls, classname):
            if classname not in cls.instances:
                cls.instances[classname] = {}
            return cls.instances[classname]

    def wrapped_calc(group_info: sgtbx.space_group_info, bravais_t: str):
        cache = MultiClassCache(bravais_t)
        info_str = group_info.type().lookup_symbol()
        try:
            result = cache[info_str]
        except KeyError:
            result = fn(group_info, bravais_t)
            cache[info_str] = result
        return result

    return wrapped_calc


@groups_cache
def calc_acentric_subgroups(
    lattice_group_info: sgtbx.space_group_info, target_bravais_str: str
):
    # Get list of sub-spacegroups
    subgrs = subgroups.subgroups(lattice_group_info).groups_parent_setting()

    # Order sub-groups
    sort_values = flex.double()
    for group in subgrs:
        order_z = group.order_z()
        space_group_number = sgtbx.space_group_type(group, False).number()
        assert 1 <= space_group_number <= 230
        sort_values.append(order_z * 1000 + space_group_number)
    perm = flex.sort_permutation(sort_values, reverse=True)
    acentric_subgroups = []
    acentric_supergroups = []
    cb_ops = []
    for i_perm in perm:
        acentric_subgroup = subgrs[i_perm]
        cb_op = acentric_subgroup.info().type().cb_op()
        ref_acentric_subgroup = acentric_subgroup.change_basis(cb_op)
        # Ignore unwanted groups
        if str(bravais_lattice(group=ref_acentric_subgroup)) != target_bravais_str:
            continue
        acentric_subgroups.append(acentric_subgroup)
        cb_ops.append(cb_op)
        acentric_supergroups.append(metric_supergroup(acentric_subgroup))
    return acentric_subgroups, acentric_supergroups, cb_ops


def find_matching_symmetry(
    unit_cell,
    target_space_group,
    max_delta=5,
    best_monoclinic_beta=True,
    target_bravais_str=None,
):
    cs = crystal.symmetry(unit_cell=unit_cell, space_group=sgtbx.space_group())
    if target_bravais_str is None:
        target_bravais_str = str(
            bravais_lattice(group=target_space_group.info().reference_setting().group())
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

    lattice_group_info = lattice_group.info()
    acentric_subgroups, acentric_supergroups, cb_ops = calc_acentric_subgroups(
        lattice_group_info,
        target_bravais_str,
    )

    for acentric_subgroup, acentric_supergroup, cb_op_minimum_ref in zip(
        acentric_subgroups, acentric_supergroups, cb_ops
    ):

        # Make symmetry object: unit-cell + space-group
        # The unit cell is potentially modified to be exactly compatible
        # with the space group symmetry.
        subsym = crystal.symmetry(
            unit_cell=minimum_symmetry.unit_cell(),
            space_group=acentric_subgroup,
            assert_is_compatible_unit_cell=False,
        )
        # Convert subgroup to reference setting
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


class SymmetryHandler:
    def __init__(self, unit_cell=None, space_group=None, max_delta=5):

        self._max_delta = max_delta
        self.target_symmetry_primitive = None
        self.target_symmetry_reference_setting = None
        self.cb_op_inp_ref = None
        self.cb_op_inp_best = None

        target_space_group = space_group
        if target_space_group is not None:
            target_space_group = target_space_group.build_derived_patterson_group()

        if unit_cell is not None:

            assert (
                space_group
            ), "space_group must be provided in combination with unit_cell"

            self.target_symmetry_inp = crystal.symmetry(
                unit_cell=unit_cell, space_group=target_space_group
            )
            if str(bravais_lattice(group=target_space_group)) == "aP":
                self.cb_op_inp_ref = (
                    self.target_symmetry_inp.change_of_basis_op_to_minimum_cell()
                )
            else:
                self.cb_op_inp_ref = (
                    self.target_symmetry_inp.change_of_basis_op_to_reference_setting()
                )
            self.target_symmetry_reference_setting = (
                self.target_symmetry_inp.change_basis(self.cb_op_inp_ref)
            )
            self.cb_op_inp_best = (
                self.target_symmetry_reference_setting.change_of_basis_op_to_best_cell()
                * self.cb_op_inp_ref
            )

        elif target_space_group is not None:
            self.target_symmetry_inp = crystal.symmetry(space_group=target_space_group)
            self.cb_op_inp_ref = (
                target_space_group.info().change_of_basis_op_to_reference_setting()
            )
            self.target_symmetry_reference_setting = crystal.symmetry(
                space_group=target_space_group.change_basis(self.cb_op_inp_ref)
            )

        cb_op_reference_to_primitive = (
            self.target_symmetry_reference_setting.change_of_basis_op_to_primitive_setting()
        )
        if unit_cell:
            self.target_symmetry_primitive = (
                self.target_symmetry_reference_setting.change_basis(
                    cb_op_reference_to_primitive
                )
            )
        else:
            self.target_symmetry_primitive = crystal.symmetry(
                space_group=self.target_symmetry_reference_setting.space_group().change_basis(
                    cb_op_reference_to_primitive
                )
            )
        self.cb_op_ref_inp = self.cb_op_inp_ref.inverse()
        self.cb_op_primitive_inp = (
            self.cb_op_ref_inp * cb_op_reference_to_primitive.inverse()
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
        logger.debug("cb_op primitive->input: %s", self.cb_op_primitive_inp)

    def apply_symmetry(self, crystal_model):
        """Apply symmetry constraints to a crystal model.

        Returns the crystal model (with symmetry constraints applied) in the same
        setting as provided as input. The cb_op returned by the method is
        that necessary to transform that model to the user-provided
        target symmetry.

        Args:
            crystal_model (dxtbx.model.Crystal): The input crystal model to which to
              apply symmetry constraints.

        Returns: (dxtbx.model.Crystal, cctbx.sgtbx.change_of_basis_op):
        The crystal model with symmetry constraints applied, and the change_of_basis_op
        that transforms the returned model to the user-specified target symmetry.
        """
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
            if bravais_lattice(group=target_sg_ref) != item["bravais"]:
                continue
            if item["max_angular_difference"] < best_angular_difference:
                best_angular_difference = item["max_angular_difference"]
                best_subgroup = item

        if best_subgroup is None:
            return None, None

        cb_op_inp_best = best_subgroup["cb_op_inp_best"]
        best_subsym = best_subgroup["best_subsym"]
        ref_subsym = best_subgroup["ref_subsym"]
        cb_op_ref_best = ref_subsym.change_of_basis_op_to_best_cell()
        cb_op_best_ref = cb_op_ref_best.inverse()
        cb_op_inp_ref = cb_op_best_ref * cb_op_inp_best
        cb_op_ref_inp = cb_op_inp_ref.inverse()

        orient = crystal_orientation(A, True)
        orient_ref = orient.change_basis(
            scitbx.matrix.sqr((cb_op_inp_ref).c().as_double_array()[0:9]).transpose()
        )
        constrain_orient = orient_ref.constrain(best_subgroup["system"])
        direct_matrix = constrain_orient.direct_matrix()

        a = scitbx.matrix.col(direct_matrix[:3])
        b = scitbx.matrix.col(direct_matrix[3:6])
        c = scitbx.matrix.col(direct_matrix[6:9])
        model = Crystal(a, b, c, space_group=target_sg_ref)
        assert target_sg_ref.is_compatible_unit_cell(model.get_unit_cell())

        model = model.change_basis(cb_op_ref_inp)

        if self.cb_op_inp_best is not None:
            # Then the unit cell has been provided: this is the cb_op to map to the
            # user-provided input unit cell
            return model, self.cb_op_inp_best.inverse() * cb_op_inp_best
        if not self.cb_op_ref_inp.is_identity_op():
            if self.target_symmetry_inp.space_group() == best_subsym.space_group():
                # Handle where e.g. the user has requested I2 instead of the reference C2
                return model, cb_op_inp_best
            # The user has specified a setting that is not the reference setting
            return model, self.cb_op_ref_inp * cb_op_inp_ref
        # Default to reference setting
        # This change of basis op will ensure that we get the best beta angle without
        # changing the centring (e.g. from C2 to I2)
        cb_op_ref_best = ref_subsym.change_of_basis_op_to_best_cell(
            best_monoclinic_beta=False
        )
        return model, cb_op_ref_best * cb_op_inp_ref
