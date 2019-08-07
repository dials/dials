#!/usr/bin/env python
# coding: utf-8
"""
Command line script to assess the space group symmetry (for MX datasets)

Currently only tests axial absences to determine which of the 65 MX space
groups is most likely.
"""
from __future__ import absolute_import, division, print_function
import logging
import sys
from dials.util import log, show_mail_on_error, Sorry
from dials.util.options import OptionParser, flatten_reflections, flatten_experiments
from dials.util.version import dials_version
from dials.util.filter_reflections import filter_reflection_table
from dials.array_family import flex
from dials.algorithms.scaling.Ih_table import (
    _reflection_table_to_iobs,
    map_indices_to_asu,
)
from dials.algorithms.symmetry.absences.laue_groups_info import (
    laue_groups,
    score_screw_axes,
    score_space_groups,
)
from dials.algorithms.symmetry.absences.screw_axes import ScrewAxisObserver
from cctbx import sgtbx
from dxtbx.model.experiment_list import ExperimentListDumper
from libtbx import phil
from libtbx.table_utils import simple_table

help_message = """Program to assess space group symmetry (for MX datasets)."""

logger = logging.getLogger("dials.space_group")
phil_scope = phil.parse(
    """
    significance_level = *0.95 0.975 0.99
        .type = choice
        .help = "Signficance to use when testing whether axial reflections are "
                "different to zero (absences and reflections in reflecting condition)"
    d_min = None
        .type = float
        .help = "Resolution limit to apply to the data."
    output {
        log = dials.systematic_absences.log
            .type = str
            .help = "The log filename"
        experiments = symmetrized.expt
            .type = str
            .help = "The name of the output file with space group updated."
        html = "absences.html"
            .type = str
            .help = "Filename for html report."
    }
"""
)


def run_sys_abs_checks(experiments, reflections, d_min=None, significance_level=0.95):
    """Check for systematic absences in the data for the laue group.

    Select the good data, merge, test screw axes and score possible space
    groups. The crystals are updated with the most likely space group.
    """

    if (
        "inverse_scale_factor" in reflections[0]
        and "intensity.scale.value" in reflections[0]
    ):
        logger.info("Attempting to perform absence checks on scaled data")
        reflections = filter_reflection_table(
            reflections[0], intensity_choice=["scale"], d_min=d_min
        )
        reflections["intensity"] = reflections["intensity.scale.value"]
        reflections["variance"] = reflections["intensity.scale.variance"]
    else:
        logger.info(
            "Attempting to perform absence checks on unscaled profile-integrated data"
        )
        reflections = filter_reflection_table(
            reflections[0], intensity_choice=["profile"], d_min=d_min
        )
        reflections["intensity"] = reflections["intensity.prf.value"]
        reflections["variance"] = reflections["intensity.prf.variance"]

    # now merge
    space_group = experiments[0].crystal.get_space_group()
    reflections["asu_miller_index"] = map_indices_to_asu(
        reflections["miller_index"], space_group
    )
    reflections["inverse_scale_factor"] = flex.double(reflections.size(), 1.0)
    merged = (
        _reflection_table_to_iobs(
            reflections, experiments[0].crystal.get_unit_cell(), space_group
        )
        .merge_equivalents(use_internal_variance=False)
        .array()
    )
    merged_reflections = flex.reflection_table()
    merged_reflections["intensity"] = merged.data()
    merged_reflections["variance"] = merged.sigmas() ** 2
    merged_reflections["miller_index"] = merged.indices()

    # Get the laue class from the space group.
    laue_group = str(space_group.build_derived_patterson_group().info())
    logger.info("Laue group: %s", laue_group)
    if laue_group not in laue_groups:
        logger.info("No absences to check for this laue group")
        return

    # Score the screw axes.
    screw_axes, screw_axis_scores = score_screw_axes(
        laue_groups[laue_group], merged_reflections, significance_level
    )

    logger.info(
        simple_table(
            [
                [
                    a.name,
                    "%.3f" % score,
                    str(a.n_refl_used[0]),
                    str(a.n_refl_used[1]),
                    "%.3f" % a.mean_I,
                    "%.3f" % a.mean_I_abs,
                    "%.3f" % a.mean_I_sigma,
                    "%.3f" % a.mean_I_sigma_abs,
                ]
                for a, score in zip(screw_axes, screw_axis_scores)
            ],
            column_headers=[
                "Screw axis",
                "Score",
                "No. present",
                "No. absent",
                "<I> present",
                "<I> absent",
                "<I/sig> present",
                "<I/sig> absent",
            ],
        ).format()
    )

    # Score the space groups from the screw axis scores.
    space_groups, scores = score_space_groups(
        screw_axis_scores, laue_groups[laue_group]
    )

    logger.info(
        simple_table(
            [[sg, "%.4f" % score] for sg, score in zip(space_groups, scores)],
            column_headers=["Space group", "score"],
        ).format()
    )

    # Find the best space group and update the experiments.
    best_sg = space_groups[scores.index(max(scores))]
    logger.info("Recommended space group: %s", best_sg)
    if "enantiomorphic pairs" in laue_groups[laue_group]:
        if best_sg in laue_groups[laue_group]["enantiomorphic pairs"]:
            logger.info(
                "Space group with equivalent score (enantiomorphic pair): %s",
                laue_groups[laue_group]["enantiomorphic pairs"][best_sg],
            )

    new_sg = sgtbx.space_group_info(symbol=best_sg).group()
    for experiment in experiments:
        experiment.crystal.set_space_group(new_sg)


def run(args=None):
    """Run the script from the command-line."""
    usage = """Usage: dials.space_group scaled.refl scaled.expt [options]"""

    parser = OptionParser(
        usage=usage,
        read_experiments=True,
        read_reflections=True,
        phil=phil_scope,
        check_format=False,
        epilog=help_message,
    )
    params, _ = parser.parse_args(args=args, show_diff_phil=False)

    if not params.input.experiments or not params.input.reflections:
        parser.print_help()
        sys.exit()

    reflections = flatten_reflections(params.input.reflections)
    experiments = flatten_experiments(params.input.experiments)

    log.config(verbosity=1, info=params.output.log)
    logger.info(dials_version())

    diff_phil = parser.diff_phil.as_str()
    if diff_phil != "":
        logger.info("The following parameters have been modified:\n")
        logger.info(diff_phil)

    ### Assert that all data have been scaled with dials - should only be
    # able to input one reflection table and experimentlist that are
    # matching and scaled together.
    if not len(reflections) == 1:
        raise Sorry("Only one reflection table can be given as input.")

    if (not "intensity.scale.value" in reflections[0]) and (
        not "intensity.prf.value" in reflections[0]
    ):
        raise Sorry(
            "Unable to find integrated or scaled reflections in the reflection table."
        )

    try:
        run_sys_abs_checks(
            experiments, reflections, params.d_min, float(params.significance_level)
        )
    except ValueError as e:
        raise Sorry(e)

    if params.output.html:
        ScrewAxisObserver().generate_html_report(params.output.html)

    if params.output.experiments:
        dump = ExperimentListDumper(experiments)
        with open(params.output.experiments, "w") as outfile:
            outfile.write(dump.as_json(split=True))


if __name__ == "__main__":
    with show_mail_on_error():
        run()
