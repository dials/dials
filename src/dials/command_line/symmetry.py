from __future__ import annotations

import copy
import json
import logging
import random
import sys

import iotbx.phil
from cctbx import sgtbx
from libtbx import Auto

import dials.util
from dials.algorithms.merging.merge import prepare_merged_reflection_table
from dials.algorithms.symmetry import (
    get_subset_for_symmetry,
    prepare_datasets_for_symmetry_analysis,
    resolution_filter_from_reflections_experiments,
)
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

exclude_inconsistent_unit_cells = False
  .type = bool
  .help = "Exclude datasets with unit cells that cannot be mapped to a common"
          "minimum cell, as controlled by the absolute_angle_tolerance and"
          "relative_length_tolerance parameters. If False, an error will be"
          "raised instead."

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

  method = *direct fourier
    .type = choice
    .help = "Use fourier analysis or direct analysis of I/sigma to determine"
            "likelihood of systematic absences"

  significance_level = 0.95
    .type = float(value_min=0, value_max=1)
    .help = "Significance to use when testing whether axial reflections are "
            "different to zero (absences and reflections in reflecting condition)."

}

output {
  log = dials.symmetry.log
    .type = str
  experiments = "symmetrized.expt"
    .type = path
  reflections = "symmetrized.refl"
    .type = path
  excluded = False
    .type = bool
  excluded_prefix = "excluded"
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

        n_datasets = len(experiments)

        datasets, experiments, refls_for_sym, cb_ops = (
            prepare_datasets_for_symmetry_analysis(
                experiments,
                reflection_tables,
                params,
                outlier_rejection_after_filter=True,
                anomalous=True,
            )
        )

        if len(datasets) != n_datasets:
            raise ValueError(
                """Some datasets have no reflection after prefiltering, please check
    input data and filtering settings e.g partiality_threshold"""
            )

        # if all datasets have been through scaling, a decision about error models has
        # been made, so don't apply any further sigma correction
        apply_sigma_correction = not all(s for s in experiments.scaling_models())

        result = LaueGroupAnalysis(
            datasets,
            normalisation=params.normalisation,
            d_min=params.d_min,
            min_i_mean_over_sigma_mean=params.min_i_mean_over_sigma_mean,
            lattice_symmetry_max_delta=params.lattice_symmetry_max_delta,
            relative_length_tolerance=params.relative_length_tolerance,
            absolute_angle_tolerance=params.absolute_angle_tolerance,
            best_monoclinic_beta=params.best_monoclinic_beta,
            apply_sigma_correction=apply_sigma_correction,
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
        _, refls_for_sym = _reindex_experiments_reflections(
            experiments, refls_for_sym, best_space_group, cb_op_inp_best
        )
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
                experiments,
                joint_reflections,
                d_min,
                partiality_threshold=params.partiality_threshold,
            )
            run_systematic_absences_checks(
                experiments,
                merged_reflections,
                float(params.systematic_absences.significance_level),
                method=params.systematic_absences.method,
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
            f"Mismatched number of experiments and reflection tables found: {len(experiments)} & {len(reflections)}."
        )
    try:
        experiments, reflections = assign_unique_identifiers(experiments, reflections)
        symmetry(experiments, reflections, params=params)
    except ValueError as e:
        sys.exit(e)


if __name__ == "__main__":
    run()
