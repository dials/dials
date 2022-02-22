from __future__ import annotations

import concurrent.futures
import copy
import logging
import math

import libtbx
import libtbx.phil
import scitbx.matrix
from cctbx import crystal, miller, sgtbx
from cctbx.crystal_orientation import crystal_orientation
from cctbx.sgtbx import bravais_types
from dxtbx.model import Crystal
from libtbx.introspection import number_of_processors
from rstbx.dps_core.lepage import iotbx_converter
from rstbx.symmetry.subgroup import MetricSubgroup
from scitbx.array_family import flex

import dials.util
from dials.algorithms.indexing.refinement import refine
from dials.command_line.check_indexing_symmetry import (
    get_symop_correlation_coefficients,
)
from dials.util.log import LoggingContext

logger = logging.getLogger(__name__)


phil_scope = libtbx.phil.parse(
    """
lepage_max_delta = 5
  .type = float
nproc = Auto
  .type = int(value_min=1)
cc_n_bins = None
  .type = int(value_min=1)
  .help = "Number of resolution bins to use for calculation of correlation coefficients"

best_monoclinic_beta = True
  .type = bool
  .help = "If True, then for monoclinic centered cells, I2 will be preferred over C2 if"
          "it gives a less oblique cell (i.e. smaller beta angle)."

include scope dials.algorithms.refinement.refiner.phil_scope
""",
    process_includes=True,
)

# override default refinement parameters
phil_scope = phil_scope.fetch(
    source=libtbx.phil.parse(
        """\
refinement {
  reflections {
    reflections_per_degree=100
  }
}
"""
    )
)


def dxtbx_crystal_from_orientation(crystal_orientation, space_group):
    """Convert a cctbx crystal_orientation to a dxtbx.Crystal model.

    Args:
        crystal_orientation (cctbx.crystal_orientation):
            A cctbx crystal_orientation object
        space_group (cctbx.sgtbx.space_group): The space group.

    Returns:
        dxtbx.model.Crystal: The dxtbx crystal model.
    """
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

    def __str__(self):
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
                min_max_cc_str = f"{item.min_cc:.3f}/{item.max_cc:.3f}"
            if item.recommended:
                status = "*"
            else:
                status = ""
            table_data.append(
                [
                    "%1s%7d" % (status, item.setting_number),
                    f"{item['max_angular_difference']:6.4f}",
                    f"{item.rmsd:5.3f}",
                    min_max_cc_str,
                    "%d" % item.Nmatches,
                    f"{item['bravais']}",
                    "%6.2f %6.2f %6.2f %6.2f %6.2f %6.2f" % P,
                    f"{uc.volume():.0f}",
                    f"{item['cb_op_inp_best'].as_abc()}",
                ]
            )

        output = dials.util.tabulate(
            table_data, headers="firstrow", colalign=("right",)
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


def lowest_symmetry_space_group_for_bravais_lattice(
    bravais_lattice: str,
) -> sgtbx.space_group:
    if bravais_lattice == "mI":
        return sgtbx.space_group_info("I2").group()
    return sgtbx.space_group_info(
        number=bravais_lattice_to_lowest_symmetry_spacegroup_number[bravais_lattice]
    ).group()


def refined_settings_from_refined_triclinic(experiments, reflections, params):
    """Generate a RefinedSettingsList from a triclinic model.

    Args:
        experiments: The experiments refined with a triclinic model
        reflections: A reflection table containing observed centroids
        params: The working PHIL parameters.

    Returns:
        RefinedSettingsList: A list of the refined settings. The highest symmetry
        setting will be first item in the list, and the triclinic setting will be last.
    """

    if params.nproc is libtbx.Auto:
        params.nproc = number_of_processors()

    if params.refinement.reflections.outlier.algorithm in ("auto", libtbx.Auto):
        if experiments[0].goniometer is None:
            params.refinement.reflections.outlier.algorithm = "sauter_poon"
        else:
            # different default to dials.refine
            # tukey is faster and more appropriate at the indexing step
            params.refinement.reflections.outlier.algorithm = "tukey"

    assert len(experiments.crystals()) == 1
    crystal = experiments.crystals()[0]

    used_reflections = copy.deepcopy(reflections)
    UC = crystal.get_unit_cell()

    refined_settings = RefinedSettingsList()
    for item in iotbx_converter(
        UC, params.lepage_max_delta, best_monoclinic_beta=params.best_monoclinic_beta
    ):
        refined_settings.append(BravaisSetting(item))

    triclinic = refined_settings.triclinic()

    # assert no transformation between indexing and bravais list
    assert str(triclinic["cb_op_inp_best"]) == "a,b,c"

    Nset = len(refined_settings)
    for j in range(Nset):
        refined_settings[j].setting_number = Nset - j

    for subgroup in refined_settings:
        bravais_lattice = str(
            bravais_types.bravais_lattice(group=subgroup["best_subsym"].space_group())
        )
        space_group = lowest_symmetry_space_group_for_bravais_lattice(bravais_lattice)
        orient = crystal_orientation(crystal.get_A(), True).change_basis(
            scitbx.matrix.sqr(
                subgroup["cb_op_inp_best"].c().as_double_array()[0:9]
            ).transpose()
        )
        constrain_orient = orient.constrain(subgroup["system"])
        subgroup["bravais"] = bravais_lattice
        subgroup.unrefined_crystal = dxtbx_crystal_from_orientation(
            constrain_orient, space_group
        )

    with concurrent.futures.ProcessPoolExecutor(max_workers=params.nproc) as pool:
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
    """Identify likely solutions using heuristics.

    Args:
        all_solutions (RefinedSettingsList): The list of refined bravais settings.

    Use a set of heuristics to identify likely solutions, by comparing refined rmsds
    in a given setting with the triclinic rmsds. Also looks at the max_angular
    difference and the correlation coefficients for the solutions. Sets the
    `recommended` attribute of each solution to `True` or `False` as appropriate.
    """
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
    higher_symmetry_miller = subgroup["cb_op_inp_best"].apply(triclinic_miller)
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
                return subgroup
            else:
                raise
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
