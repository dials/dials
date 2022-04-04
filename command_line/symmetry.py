from __future__ import annotations

import collections
import copy
import json
import logging
import math
import random
import sys

import iotbx.phil
from cctbx import sgtbx, uctbx
from cctbx.sgtbx.bravais_types import bravais_lattice
from cctbx.sgtbx.lattice_symmetry import metric_subgroups
from dxtbx.model import ExperimentList
from libtbx import Auto

import dials.util
from dials.algorithms.merging.merge import prepare_merged_reflection_table
from dials.algorithms.symmetry import resolution_filter_from_reflections_experiments
from dials.algorithms.symmetry.absences.laue_groups_info import (
    laue_groups as laue_groups_for_absence_analysis,
)
from dials.algorithms.symmetry.absences.run_absences_checks import (
    run_systematic_absences_checks,
)
from dials.algorithms.symmetry.absences.screw_axes import ScrewAxisObserver
from dials.algorithms.symmetry.laue_group import LaueGroupAnalysis
from dials.array_family import flex
from dials.command_line.reindex import reindex_experiments
from dials.util import log, tabulate
from dials.util.exclude_images import (
    exclude_image_ranges_from_scans,
    get_selection_for_valid_image_ranges,
)
from dials.util.filter_reflections import filtered_arrays_from_experiments_reflections
from dials.util.multi_dataset_handling import (
    assign_unique_identifiers,
    parse_multiple_datasets,
    update_imageset_ids,
)
from dials.util.options import ArgumentParser, reflections_and_experiments_from_files
from dials.util.version import dials_version

logger = logging.getLogger("dials.command_line.symmetry")

phil_scope = iotbx.phil.parse(
    """\
include scope dials.util.exclude_images.phil_scope
d_min = Auto
  .type = float(value_min=0)

min_i_mean_over_sigma_mean = 4
  .type = float(value_min=0)

min_cc_half = 0.6
  .type = float(value_min=0, value_max=1)

normalisation = kernel quasi ml_iso *ml_aniso
  .type = choice

lattice_group = None
  .type = space_group

seed = 230
  .type = int(value_min=0)

lattice_symmetry_max_delta = 2.0
  .type = float(value_min=0)

relative_length_tolerance = 0.05
  .type = float(value_min=0)

absolute_angle_tolerance = 2
  .type = float(value_min=0)

partiality_threshold = 0.4
  .type = float
  .help = "Use only reflections with a partiality above this threshold."

laue_group = auto
  .type = space_group
  .help = "Optionally specify the Laue group. If set to auto, then test all possible "
          "Laue groups. If set to None, then take the Laue group from the input file."

change_of_basis_op = None
  .type = str

best_monoclinic_beta = True
  .type = bool
  .help = "If True, then for monoclinic centered cells, I2 will be preferred over C2 if"
          "it gives a less oblique cell (i.e. smaller beta angle)."

systematic_absences {

  check = True
    .type = bool
    .help = "Check systematic absences for the current laue group."

  significance_level = *0.95 0.975 0.99
    .type = choice
    .help = "Signficance to use when testing whether axial reflections are "
            "different to zero (absences and reflections in reflecting condition)."

}

output {
  log = dials.symmetry.log
    .type = str
  experiments = "symmetrized.expt"
    .type = path
  reflections = "symmetrized.refl"
    .type = path
  json = dials.symmetry.json
    .type = path
  html = "dials.symmetry.html"
    .type = path
    .help = "Filename for html report."
}
""",
    process_includes=True,
)


def median_unit_cell(experiments):
    uc_params = [flex.double() for i in range(6)]
    for c in experiments.crystals():
        for i, p in enumerate(c.get_unit_cell().parameters()):
            uc_params[i].append(p)
    return uctbx.unit_cell(parameters=[flex.median(p) for p in uc_params])


def unit_cells_are_similar_to(
    experiments, unit_cell, relative_length_tolerance, absolute_angle_tolerance
):
    return all(
        expt.crystal.get_unit_cell().is_similar_to(
            unit_cell,
            relative_length_tolerance=relative_length_tolerance,
            absolute_angle_tolerance=absolute_angle_tolerance,
        )
        for expt in experiments
    )


def change_of_basis_ops_to_minimum_cell(
    experiments, max_delta, relative_length_tolerance, absolute_angle_tolerance
):
    """
    Compute change of basis ops to map experiments to the minimum cell

    Map to the minimum cell via the best cell, which appears to guarantee that the
    resulting minimum cells are consistent.

    Args:
        experiments (ExperimentList): a list of experiments.
        reflections (list): a list of reflection tables

    Returns: The experiments and reflections mapped to the minimum cell
    """

    median_cell = median_unit_cell(experiments)
    unit_cells_are_similar = unit_cells_are_similar_to(
        experiments, median_cell, relative_length_tolerance, absolute_angle_tolerance
    )
    centring_symbols = [
        bravais_lattice(group=expt.crystal.get_space_group()).centring_symbol
        for expt in experiments
    ]
    if unit_cells_are_similar and len(set(centring_symbols)) == 1:
        groups = metric_subgroups(
            experiments[0]
            .crystal.get_crystal_symmetry()
            .customized_copy(unit_cell=median_cell),
            max_delta,
            enforce_max_delta_for_generated_two_folds=True,
        )
        group = groups.result_groups[0]
        cb_op_best_to_min = group["best_subsym"].change_of_basis_op_to_minimum_cell()
        cb_ops = [cb_op_best_to_min * group["cb_op_inp_best"]] * len(experiments)
    else:
        groups = [
            metric_subgroups(
                expt.crystal.get_crystal_symmetry(),
                max_delta,
                best_monoclinic_beta=False,
                enforce_max_delta_for_generated_two_folds=True,
            )
            for expt in experiments
        ]
        counter = collections.Counter(
            g.result_groups[0]["best_subsym"].space_group() for g in groups
        )
        target_group = counter.most_common()[0][0]
        cb_ops = []
        for expt in experiments:
            groups = metric_subgroups(
                expt.crystal.get_crystal_symmetry(),
                max_delta,
                best_monoclinic_beta=False,
                enforce_max_delta_for_generated_two_folds=True,
            )
            group = None
            for g in groups.result_groups:
                if g["best_subsym"].space_group() == target_group:
                    group = g
            if group:
                cb_ops.append(group["cb_op_inp_best"])
            else:
                cb_ops.append(None)
                logger.info(
                    f"Couldn't match unit cell to target symmetry:\n"
                    f"{expt.crystal.get_crystal_symmetry()}\n"
                    f"Target symmetry: {target_group.info()}"
                )
        ref_expts = ExperimentList(
            [expt for expt, cb_op in zip(experiments, cb_ops) if cb_op]
        ).change_basis(list(filter(None, cb_ops)))
        cb_op_ref_min = (
            ref_expts[0]
            .crystal.get_crystal_symmetry()
            .customized_copy(unit_cell=median_unit_cell(ref_expts))
            .change_of_basis_op_to_minimum_cell()
        )
        cb_ops = [cb_op_ref_min * cb_op if cb_op else None for cb_op in cb_ops]
    return cb_ops


def apply_change_of_basis_ops(experiments, reflections, change_of_basis_ops):
    """
    Apply the given change of basis ops to the input experiments and reflections

    Args:
        experiments (ExperimentList): a list of experiments.
        reflections (list): a list of reflection tables
        change_of_basis_ops (list): a list of cctbx.sgtbx.change_of_basis_op

    Returns: The experiments and reflections after application of the change of basis ops
    """

    for expt, refl, cb_op_inp_min in zip(experiments, reflections, change_of_basis_ops):
        refl["miller_index"] = cb_op_inp_min.apply(refl["miller_index"])
        expt.crystal = expt.crystal.change_basis(cb_op_inp_min)
        expt.crystal.set_space_group(sgtbx.space_group())
    return experiments, reflections


def eliminate_sys_absent(experiments, reflections):
    for i, expt in enumerate(experiments):
        if expt.crystal.get_space_group().n_ltr() > 1:
            effective_group = (
                expt.crystal.get_space_group().build_derived_reflection_intensity_group(
                    anomalous_flag=True
                )
            )
            sys_absent_flags = effective_group.is_sys_absent(
                reflections[i]["miller_index"]
            )
            if sys_absent_flags.count(True):
                reflections[i] = reflections[i].select(~sys_absent_flags)
                logger.info(
                    "Eliminating %i systematic absences for experiment %s",
                    sys_absent_flags.count(True),
                    expt.identifier,
                )
    return reflections


def get_subset_for_symmetry(experiments, reflection_tables, exclude_images=None):
    """Select an image range for symmetry analysis, or just select
    the first 360 degrees of data."""
    refls_for_sym = []
    if exclude_images:
        experiments = exclude_image_ranges_from_scans(
            reflection_tables, experiments, exclude_images
        )
        for refl, exp in zip(reflection_tables, experiments):
            sel = get_selection_for_valid_image_ranges(refl, exp)
            refls_for_sym.append(refl.select(sel))
    else:
        for expt, refl in zip(experiments, reflection_tables):
            sel = get_selection_for_valid_image_ranges(refl, expt)
            if not sel.count(False):
                # Use first 360 degrees if <360 deg i.e. first measured data,
                # but only if no reflections have been explicitly excluded
                # already
                scan_end = int(math.ceil(360 / abs(expt.scan.get_oscillation()[1])))
                if scan_end < len(expt.scan):
                    sel = refl["xyzobs.px.value"].parts()[2] <= scan_end
            refls_for_sym.append(refl.select(sel))
    return refls_for_sym


def symmetry(experiments, reflection_tables, params=None):
    """
    Run symmetry analysis

    Args:
        experiments: An experiment list.
        reflection_tables: A list of reflection tables.
        params: The dials.symmetry phil scope.
    """
    result = None
    if params is None:
        params = phil_scope.extract()
    refls_for_sym = []

    if params.laue_group is Auto:
        logger.info("=" * 80)
        logger.info("")
        logger.info("Performing Laue group analysis")
        logger.info("")

        # Transform models into miller arrays
        n_datasets = len(experiments)

        # Map experiments and reflections to minimum cell
        # Eliminate reflections that are systematically absent due to centring
        # of the lattice, otherwise they would lead to non-integer miller indices
        # when reindexing to a primitive setting
        cb_ops = change_of_basis_ops_to_minimum_cell(
            experiments,
            params.lattice_symmetry_max_delta,
            params.relative_length_tolerance,
            params.absolute_angle_tolerance,
        )
        reflection_tables = eliminate_sys_absent(experiments, reflection_tables)
        experiments, reflection_tables = apply_change_of_basis_ops(
            experiments, reflection_tables, cb_ops
        )

        refls_for_sym = get_subset_for_symmetry(
            experiments, reflection_tables, params.exclude_images
        )

        datasets = filtered_arrays_from_experiments_reflections(
            experiments,
            refls_for_sym,
            outlier_rejection_after_filter=True,
            partiality_threshold=params.partiality_threshold,
        )
        if len(datasets) != n_datasets:
            raise ValueError(
                """Some datasets have no reflection after prefiltering, please check
    input data and filtering settings e.g partiality_threshold"""
            )

        datasets = [
            ma.as_anomalous_array().merge_equivalents().array() for ma in datasets
        ]
        result = LaueGroupAnalysis(
            datasets,
            normalisation=params.normalisation,
            d_min=params.d_min,
            min_i_mean_over_sigma_mean=params.min_i_mean_over_sigma_mean,
            lattice_symmetry_max_delta=params.lattice_symmetry_max_delta,
            relative_length_tolerance=params.relative_length_tolerance,
            absolute_angle_tolerance=params.absolute_angle_tolerance,
            best_monoclinic_beta=params.best_monoclinic_beta,
        )
        logger.info("")
        logger.info(result)

        if params.output.json is not None:
            d = result.as_dict()
            d["cb_op_inp_min"] = [str(cb_op) for cb_op in cb_ops]
            # Copy the "input_symmetry" to "min_cell_symmetry" as it isn't technically
            # the input symmetry to dials.symmetry
            d["min_cell_symmetry"] = d["input_symmetry"]
            del d["input_symmetry"]
            json_str = json.dumps(d, indent=2)
            with open(params.output.json, "w") as f:
                f.write(json_str)

        # Change of basis operator from input unit cell to best unit cell
        cb_op_inp_best = result.best_solution.subgroup["cb_op_inp_best"]
        # Get the best space group.
        best_subsym = result.best_solution.subgroup["best_subsym"]
        best_space_group = best_subsym.space_group().build_derived_acentric_group()
        logger.info(
            tabulate(
                [[str(best_subsym.space_group_info()), str(best_space_group.info())]],
                ["Patterson group", "Corresponding MX group"],
            )
        )
        # Reindex the input data
        experiments, reflection_tables = _reindex_experiments_reflections(
            experiments, reflection_tables, best_space_group, cb_op_inp_best
        )

    elif params.laue_group is not None:
        if params.change_of_basis_op is not None:
            cb_op = sgtbx.change_of_basis_op(params.change_of_basis_op)
        else:
            cb_op = sgtbx.change_of_basis_op()
        # Reindex the input data
        experiments, reflection_tables = _reindex_experiments_reflections(
            experiments, reflection_tables, params.laue_group.group(), cb_op
        )

    if params.systematic_absences.check:
        logger.info("=" * 80)
        logger.info("")
        logger.info("Analysing systematic absences")
        logger.info("")

        # Get the laue class from the current space group.
        space_group = experiments[0].crystal.get_space_group()
        laue_group = str(space_group.build_derived_patterson_group().info())
        logger.info("Laue group: %s", laue_group)
        if laue_group in ("I m -3", "I m m m"):
            if laue_group == "I m -3":
                logger.info(
                    """Space groups I 2 3 & I 21 3 cannot be distinguished with systematic absence
analysis, due to lattice centering.
Using space group I 2 3, space group I 21 3 is equally likely.\n"""
                )
            if laue_group == "I m m m":
                logger.info(
                    """Space groups I 2 2 2 & I 21 21 21 cannot be distinguished with systematic absence
analysis, due to lattice centering.
Using space group I 2 2 2, space group I 21 21 21 is equally likely.\n"""
                )
        elif laue_group not in laue_groups_for_absence_analysis:
            logger.info("No absences to check for this laue group\n")
        else:
            if not refls_for_sym:
                refls_for_sym = get_subset_for_symmetry(
                    experiments, reflection_tables, params.exclude_images
                )

            if (params.d_min is Auto) and (result is not None):
                d_min = result.intensities.resolution_range()[1]
            elif params.d_min is Auto:
                d_min = resolution_filter_from_reflections_experiments(
                    refls_for_sym,
                    experiments,
                    params.min_i_mean_over_sigma_mean,
                    params.min_cc_half,
                )
            else:
                d_min = params.d_min

            # combine before sys abs test - only triggers if laue_group=None and
            # multiple input files.
            if len(reflection_tables) > 1:
                joint_reflections = flex.reflection_table()
                for table in refls_for_sym:
                    joint_reflections.extend(table)
            else:
                joint_reflections = refls_for_sym[0]

            merged_reflections = prepare_merged_reflection_table(
                experiments, joint_reflections, d_min
            )
            run_systematic_absences_checks(
                experiments,
                merged_reflections,
                float(params.systematic_absences.significance_level),
            )

    logger.info(
        "Saving reindexed experiments to %s in space group %s",
        params.output.experiments,
        str(experiments[0].crystal.get_space_group().info()),
    )
    experiments.as_file(params.output.experiments)
    if params.output.reflections is not None:
        if len(reflection_tables) > 1:
            joint_reflections = flex.reflection_table()
            for table in reflection_tables:
                joint_reflections.extend(table)
        else:
            joint_reflections = reflection_tables[0]
        logger.info(
            "Saving %s reindexed reflections to %s",
            len(joint_reflections),
            params.output.reflections,
        )
        joint_reflections.as_file(params.output.reflections)

    if params.output.html and params.systematic_absences.check:
        ScrewAxisObserver().generate_html_report(params.output.html)


def _reindex_experiments_reflections(experiments, reflections, space_group, cb_op):
    """Reindex the input data."""
    reindexed_experiments = reindex_experiments(
        experiments, cb_op, space_group=space_group
    )
    reindexed_reflections = flex.reflection_table()
    reflections = update_imageset_ids(experiments, reflections)
    for i in range(len(reindexed_experiments)):
        reindexed_refl = copy.deepcopy(reflections[i])
        reindexed_refl["miller_index"] = cb_op.apply(reindexed_refl["miller_index"])
        reindexed_reflections.extend(reindexed_refl)
    return reindexed_experiments, [reindexed_reflections]


help_message = """
This program implements the methods of
`POINTLESS <http://www.ccp4.ac.uk/html/pointless.html>`_ (
`Evans, P. (2006). Acta Cryst. D62, 72-82. <https://doi.org/10.1107/S0907444905036693>`_ and
`Evans, P. R. (2011). Acta Cryst. D67, 282-292. <https://doi.org/10.1107/S090744491003982X>`_)
for scoring and determination of Laue group symmetry.

The program takes as input a set of one or more integrated experiments and
reflections.

Examples::

  dials.symmetry models.expt observations.refl
"""


@dials.util.show_mail_handle_errors()
def run(args=None):
    """Run symmetry analysis from the command-line."""
    usage = "dials.symmetry [options] models.expt observations.refl"

    parser = ArgumentParser(
        usage=usage,
        phil=phil_scope,
        read_reflections=True,
        read_experiments=True,
        check_format=False,
        epilog=help_message,
    )

    params, options, args = parser.parse_args(
        args=args, show_diff_phil=False, return_unhandled=True
    )

    # Configure the logging
    log.config(verbosity=options.verbose, logfile=params.output.log)

    logger.info(dials_version())

    # Log the diff phil
    diff_phil = parser.diff_phil.as_str()
    if diff_phil != "":
        logger.info("The following parameters have been modified:\n")
        logger.info(diff_phil)

    if params.seed is not None:
        flex.set_random_seed(params.seed)
        random.seed(params.seed)

    if not params.input.experiments or not params.input.reflections:
        parser.print_help()
        sys.exit()

    reflections, experiments = reflections_and_experiments_from_files(
        params.input.reflections, params.input.experiments
    )

    reflections = parse_multiple_datasets(reflections)

    if len(experiments) != len(reflections):
        sys.exit(
            "Mismatched number of experiments and reflection tables found: %s & %s."
            % (len(experiments), len(reflections))
        )
    try:
        experiments, reflections = assign_unique_identifiers(experiments, reflections)
        symmetry(experiments, reflections, params=params)
    except ValueError as e:
        sys.exit(e)


if __name__ == "__main__":
    run()
