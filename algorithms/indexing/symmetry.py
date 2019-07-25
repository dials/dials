from __future__ import absolute_import, division, print_function

import copy
import logging
import math

import scitbx.matrix
from cctbx import crystal, sgtbx
from cctbx.crystal_orientation import crystal_orientation
from cctbx.sgtbx import bravais_types, change_of_basis_op, subgroups
from cctbx.sgtbx import lattice_symmetry
from cctbx.sgtbx.bravais_types import bravais_lattice
from dials.algorithms.indexing import DialsIndexError
from dxtbx.model import Crystal
from libtbx import easy_mp
from rstbx.dps_core.lepage import iotbx_converter
from rstbx.symmetry.subgroup import MetricSubgroup
from scitbx.array_family import flex
from six.moves import StringIO

logger = logging.getLogger(__name__)


def dials_crystal_from_orientation(crystal_orientation, space_group):
    dm = crystal_orientation.direct_matrix()
    AA = scitbx.matrix.col((dm[0], dm[1], dm[2]))
    BB = scitbx.matrix.col((dm[3], dm[4], dm[5]))
    CC = scitbx.matrix.col((dm[6], dm[7], dm[8]))

    cryst = Crystal(
        real_space_a=AA, real_space_b=BB, real_space_c=CC, space_group=space_group
    )
    return cryst


class BravaisSetting(MetricSubgroup):  # inherits from dictionary
    def __init__(self, other):
        self.update(other)


class RefinedSettingsList(list):
    def supergroup(self):
        return self[0]

    def triclinic(self):
        return self[-1]

    def as_dict(self):
        result = {}

        for item in self:
            uc = item.refined_crystal.get_unit_cell()
            result[item.setting_number] = {
                "max_angular_difference": item["max_angular_difference"],
                "rmsd": item.rmsd,
                "nspots": item.Nmatches,
                "bravais": item["bravais"],
                "unit_cell": uc.parameters(),
                "cb_op": item["cb_op_inp_best"].as_abc(),
                "max_cc": item.max_cc,
                "min_cc": item.min_cc,
                "correlation_coefficients": list(item.correlation_coefficients),
                "cc_nrefs": list(item.cc_nrefs),
                "recommended": item.recommended,
            }

        return result

    def labelit_printout(self):
        from libtbx import table_utils

        table_data = [
            [
                "Solution",
                "Metric fit",
                "rmsd",
                "min/max cc",
                "#spots",
                "lattice",
                "unit_cell",
                "volume",
                "cb_op",
            ]
        ]
        for item in self:
            uc = item.refined_crystal.get_unit_cell()
            P = uc.parameters()
            min_max_cc_str = "-/-"
            if item.min_cc is not None and item.max_cc is not None:
                min_max_cc_str = "%.3f/%.3f" % (item.min_cc, item.max_cc)
            if item.recommended:
                status = "*"
            else:
                status = ""
            table_data.append(
                [
                    "%1s%7d" % (status, item.setting_number),
                    "%(max_angular_difference)6.4f" % item,
                    "%5.3f" % item.rmsd,
                    min_max_cc_str,
                    "%d" % item.Nmatches,
                    "%(bravais)s" % item,
                    "%6.2f %6.2f %6.2f %6.2f %6.2f %6.2f" % P,
                    "%.0f" % uc.volume(),
                    "%s" % item["cb_op_inp_best"].as_abc(),
                ]
            )

        output = table_utils.format(
            table_data, has_header=1, justify="right", delim=" "
        )
        output = output + "\n* = recommended solution\n"
        return output


# Mapping of Bravais lattice type to corresponding lowest possible symmetry
bravais_lattice_to_lowest_symmetry_spacegroup_number = {
    "aP": 1,
    "mP": 3,
    "mC": 5,
    "oP": 16,
    "oC": 20,
    "oF": 22,
    "oI": 23,
    "tP": 75,
    "tI": 79,
    "hP": 143,
    "hR": 146,
    "cP": 195,
    "cF": 196,
    "cI": 197,
}


def refined_settings_factory_from_refined_triclinic(
    params,
    experiments,
    reflections,
    i_setting=None,
    lepage_max_delta=5.0,
    nproc=1,
    refiner_verbosity=0,
):

    assert len(experiments.crystals()) == 1
    crystal = experiments.crystals()[0]

    used_reflections = copy.deepcopy(reflections)
    UC = crystal.get_unit_cell()

    Lfat = RefinedSettingsList()
    for item in iotbx_converter(UC, lepage_max_delta):
        Lfat.append(BravaisSetting(item))

    triclinic = Lfat.triclinic()

    # assert no transformation between indexing and bravais list
    assert str(triclinic["cb_op_inp_best"]) == "a,b,c"

    Nset = len(Lfat)
    for j in range(Nset):
        Lfat[j].setting_number = Nset - j

    for j in range(Nset):
        cb_op = Lfat[j]["cb_op_inp_best"].c().as_double_array()[0:9]
        orient = crystal_orientation(crystal.get_A(), True)
        orient_best = orient.change_basis(scitbx.matrix.sqr(cb_op).transpose())
        constrain_orient = orient_best.constrain(Lfat[j]["system"])
        bravais = Lfat[j]["bravais"]
        cb_op_best_ref = Lfat[j][
            "best_subsym"
        ].change_of_basis_op_to_reference_setting()
        space_group = sgtbx.space_group_info(
            number=bravais_lattice_to_lowest_symmetry_spacegroup_number[bravais]
        ).group()
        space_group = space_group.change_basis(cb_op_best_ref.inverse())
        bravais = str(bravais_types.bravais_lattice(group=space_group))
        Lfat[j]["bravais"] = bravais
        Lfat[j].unrefined_crystal = dials_crystal_from_orientation(
            constrain_orient, space_group
        )

    args = []
    for subgroup in Lfat:
        args.append(
            (params, subgroup, used_reflections, experiments, refiner_verbosity)
        )

    results = easy_mp.parallel_map(
        func=refine_subgroup,
        iterable=args,
        processes=nproc,
        method="multiprocessing",
        preserve_order=True,
        asynchronous=True,
        preserve_exception_message=True,
    )

    for i, result in enumerate(results):
        Lfat[i] = result
    identify_likely_solutions(Lfat)
    return Lfat


def identify_likely_solutions(all_solutions):
    p1_solution = all_solutions[-1]
    assert p1_solution.setting_number == 1, p1_solution.setting_number
    rmsd_p1 = p1_solution.rmsd

    for solution in all_solutions:
        solution.recommended = False
        if solution["max_angular_difference"] < 0.5:
            if (
                solution.min_cc is None or solution.min_cc < 0.5
            ) and solution.rmsd > 1.5 * rmsd_p1:
                continue
        elif solution.min_cc < 0.7 and solution.rmsd > 2.0 * rmsd_p1:
            continue
        elif solution.rmsd > 3 * rmsd_p1:
            continue
        solution.recommended = True


def refine_subgroup(args):
    assert len(args) == 5
    from dials.command_line.check_indexing_symmetry import (
        get_symop_correlation_coefficients,
        normalise_intensities,
    )

    params, subgroup, used_reflections, experiments, refiner_verbosity = args

    used_reflections = copy.deepcopy(used_reflections)
    triclinic_miller = used_reflections["miller_index"]
    cb_op = subgroup["cb_op_inp_best"]
    higher_symmetry_miller = cb_op.apply(triclinic_miller)
    used_reflections["miller_index"] = higher_symmetry_miller
    unrefined_crystal = copy.deepcopy(subgroup.unrefined_crystal)
    for expt in experiments:
        expt.crystal = unrefined_crystal

    from dials.algorithms.indexing.refinement import refine

    subgroup.max_cc = None
    subgroup.min_cc = None
    subgroup.correlation_coefficients = []
    subgroup.cc_nrefs = []
    try:
        logger = logging.getLogger()
        level = logger.getEffectiveLevel()
        logger.setLevel(logging.ERROR)
        iqr_multiplier = params.refinement.reflections.outlier.tukey.iqr_multiplier
        params.refinement.reflections.outlier.tukey.iqr_multiplier = 2 * iqr_multiplier
        refinery, refined, outliers = refine(
            params, used_reflections, experiments, verbosity=refiner_verbosity
        )
        params.refinement.reflections.outlier.tukey.iqr_multiplier = iqr_multiplier
        refinery, refined, outliers = refine(
            params,
            used_reflections,
            refinery.get_experiments(),
            verbosity=refiner_verbosity,
        )
    except RuntimeError as e:
        if (
            str(e) == "scitbx Error: g0 - astry*astry -astrz*astrz <= 0."
            or str(e) == "scitbx Error: g1-bstrz*bstrz <= 0."
        ):
            subgroup.refined_experiments = None
            subgroup.rmsd = None
            subgroup.Nmatches = None
        else:
            raise
    else:
        dall = refinery.rmsds()
        dx = dall[0]
        dy = dall[1]
        subgroup.rmsd = math.sqrt(dx * dx + dy * dy)
        subgroup.Nmatches = len(refinery.get_matches())
        subgroup.refined_experiments = refinery.get_experiments()
        assert len(subgroup.refined_experiments.crystals()) == 1
        subgroup.refined_crystal = subgroup.refined_experiments.crystals()[0]
        cs = crystal.symmetry(
            unit_cell=subgroup.refined_crystal.get_unit_cell(),
            space_group=subgroup.refined_crystal.get_space_group(),
        )
        if "intensity.sum.value" in used_reflections:
            # remove refl with -ve variance
            sel = used_reflections["intensity.sum.variance"] > 0
            good_reflections = used_reflections.select(sel)
            from cctbx import miller

            ms = miller.set(cs, good_reflections["miller_index"])
            ms = ms.array(
                good_reflections["intensity.sum.value"]
                / flex.sqrt(good_reflections["intensity.sum.variance"])
            )
            if params.normalise:
                if params.normalise_bins:
                    ms = normalise_intensities(ms, n_bins=params.normalise_bins)
                else:
                    ms = normalise_intensities(ms)
            if params.cc_n_bins is not None:
                ms.setup_binner(n_bins=params.cc_n_bins)
            ccs, nrefs = get_symop_correlation_coefficients(
                ms, use_binning=(params.cc_n_bins is not None)
            )
            subgroup.correlation_coefficients = ccs
            subgroup.cc_nrefs = nrefs
            ccs = ccs.select(nrefs > 10)
            if len(ccs) > 1:
                subgroup.max_cc = flex.max(ccs[1:])
                subgroup.min_cc = flex.min(ccs[1:])
    finally:
        logger.setLevel(level)
    return subgroup


find_max_delta = sgtbx.lattice_symmetry_find_max_delta


def metric_supergroup(group):
    return (
        sgtbx.space_group_info(group=group)
        .type()
        .expand_addl_generators_of_euclidean_normalizer(True, True)
        .build_derived_acentric_group()
    )


def find_matching_symmetry(unit_cell, target_space_group, max_delta=5):
    cs = crystal.symmetry(unit_cell=unit_cell, space_group=sgtbx.space_group())
    target_bravais_t = bravais_types.bravais_lattice(
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
        bravais_t = bravais_types.bravais_lattice(group=ref_subsym.space_group())
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

        max_angular_difference = find_max_delta(
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

            if self.target_symmetry_reference_setting is not None:
                debug_handle = StringIO()
                self.target_symmetry_reference_setting.show_summary(f=debug_handle)
                logger.debug(
                    "Target symmetry (reference setting):\n" + debug_handle.getvalue()
                )
            if self.target_symmetry_primitive is not None:
                debug_handle = StringIO()
                self.target_symmetry_primitive.show_summary(f=debug_handle)
                logger.debug(
                    "Target symmetry (primitive cell):\n" + debug_handle.getvalue()
                )
            logger.debug(
                "cb_op reference->primitive: " + str(self.cb_op_reference_to_primitive)
            )
            logger.debug("cb_op primitive->input: " + str(self.cb_op_primitive_inp))

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
