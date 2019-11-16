from __future__ import absolute_import, division, print_function

import concurrent.futures
import copy
import logging
import math

import scitbx.matrix
from cctbx import crystal, miller, sgtbx
from cctbx.crystal_orientation import crystal_orientation
from cctbx.sgtbx import bravais_types
from rstbx.dps_core.lepage import iotbx_converter
from rstbx.symmetry.subgroup import MetricSubgroup
from scitbx.array_family import flex

from dxtbx.model import Crystal

import dials.util
from dials.algorithms.indexing.refinement import refine
from dials.command_line.check_indexing_symmetry import (
    get_symop_correlation_coefficients,
)
from dials.util.log import LoggingContext

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

    def as_str(self):
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

        output = dials.util.tabulate(table_data, headers="firstrow")
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


def refined_settings_from_refined_triclinic(
    params, experiments, reflections, lepage_max_delta=5.0, nproc=1
):

    assert len(experiments.crystals()) == 1
    crystal = experiments.crystals()[0]

    used_reflections = copy.deepcopy(reflections)
    UC = crystal.get_unit_cell()

    refined_settings = RefinedSettingsList()
    for item in iotbx_converter(UC, lepage_max_delta):
        refined_settings.append(BravaisSetting(item))

    triclinic = refined_settings.triclinic()

    # assert no transformation between indexing and bravais list
    assert str(triclinic["cb_op_inp_best"]) == "a,b,c"

    Nset = len(refined_settings)
    for j in range(Nset):
        refined_settings[j].setting_number = Nset - j

    for j in range(Nset):
        cb_op = refined_settings[j]["cb_op_inp_best"].c().as_double_array()[0:9]
        orient = crystal_orientation(crystal.get_A(), True)
        orient_best = orient.change_basis(scitbx.matrix.sqr(cb_op).transpose())
        constrain_orient = orient_best.constrain(refined_settings[j]["system"])
        bravais = refined_settings[j]["bravais"]
        cb_op_best_ref = refined_settings[j][
            "best_subsym"
        ].change_of_basis_op_to_reference_setting()
        space_group = sgtbx.space_group_info(
            number=bravais_lattice_to_lowest_symmetry_spacegroup_number[bravais]
        ).group()
        space_group = space_group.change_basis(cb_op_best_ref.inverse())
        bravais = str(bravais_types.bravais_lattice(group=space_group))
        refined_settings[j]["bravais"] = bravais
        refined_settings[j].unrefined_crystal = dials_crystal_from_orientation(
            constrain_orient, space_group
        )

    with concurrent.futures.ProcessPoolExecutor(max_workers=nproc) as pool:
        for i, result in enumerate(
            pool.map(
                refine_subgroup,
                (
                    (params, subgroup, used_reflections, experiments)
                    for subgroup in refined_settings
                ),
            )
        ):
            refined_settings[i] = result

    identify_likely_solutions(refined_settings)
    return refined_settings


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
        elif (
            solution.min_cc is None or solution.min_cc < 0.7
        ) and solution.rmsd > 2.0 * rmsd_p1:
            continue
        elif solution.rmsd > 3 * rmsd_p1:
            continue
        solution.recommended = True


def refine_subgroup(args):
    assert len(args) == 4
    params, subgroup, used_reflections, experiments = args

    used_reflections = copy.deepcopy(used_reflections)
    triclinic_miller = used_reflections["miller_index"]
    cb_op = subgroup["cb_op_inp_best"]
    higher_symmetry_miller = cb_op.apply(triclinic_miller)
    used_reflections["miller_index"] = higher_symmetry_miller
    unrefined_crystal = copy.deepcopy(subgroup.unrefined_crystal)
    for expt in experiments:
        expt.crystal = unrefined_crystal

    subgroup.max_cc = None
    subgroup.min_cc = None
    subgroup.correlation_coefficients = []
    subgroup.cc_nrefs = []

    with LoggingContext("dials.algorithms.refinement", level=logging.ERROR):
        try:
            outlier_algorithm = params.refinement.reflections.outlier.algorithm
            sel = used_reflections.get_flags(used_reflections.flags.used_in_refinement)
            if sel.all_eq(False):
                # Soft outlier rejection if no used_in_refinement flag is set
                params.refinement.reflections.outlier.algorithm = "tukey"
                iqr_multiplier = (
                    params.refinement.reflections.outlier.tukey.iqr_multiplier
                )
                params.refinement.reflections.outlier.tukey.iqr_multiplier = (
                    2 * iqr_multiplier
                )
                sel = ~sel
            else:
                # Remove reflections not previously used in refinement
                params.refinement.reflections.outlier.algorithm = "null"
            refinery, refined, outliers = refine(
                params, used_reflections.select(sel), experiments
            )
            params.refinement.reflections.outlier.algorithm = outlier_algorithm
            refinery, refined, outliers = refine(
                params, used_reflections, refinery.get_experiments()
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

                ms = miller.set(cs, good_reflections["miller_index"])
                ms = ms.array(
                    good_reflections["intensity.sum.value"]
                    / flex.sqrt(good_reflections["intensity.sum.variance"])
                )
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
    return subgroup
