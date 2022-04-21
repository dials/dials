"""Definition of systematic absences check algorithm."""
from __future__ import annotations

import logging

from cctbx import sgtbx

from dials.algorithms.symmetry.absences.laue_groups_info import (
    laue_groups,
    score_screw_axes,
    score_space_groups,
)
from dials.util import tabulate

logger = logging.getLogger("dials")


def run_systematic_absences_checks(
    experiments, merged_reflections, significance_level=0.95
):
    """Check for systematic absences in the data for the laue group.

    Using a reflection table containing merged data, test screw axes and score
    possible space groups. The crystals are updated with the most likely space
    group.
    """
    # Get the laue class from the space group.
    space_group = experiments[0].crystal.get_space_group()
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
        "%s",
        tabulate(
            [
                [
                    a.name,
                    f"{score:.3f}",
                    str(a.n_refl_used[0]),
                    str(a.n_refl_used[1]),
                    f"{a.mean_I:.3f}",
                    f"{a.mean_I_abs:.3f}",
                    f"{a.mean_I_sigma:.3f}",
                    f"{a.mean_I_sigma_abs:.3f}",
                ]
                for a, score in zip(screw_axes, screw_axis_scores)
            ],
            [
                "Screw axis",
                "Score",
                "No. present",
                "No. absent",
                "<I> present",
                "<I> absent",
                "<I/sig> present",
                "<I/sig> absent",
            ],
        ),
    )

    # Score the space groups from the screw axis scores.
    space_groups, scores = score_space_groups(
        screw_axis_scores, laue_groups[laue_group]
    )

    logger.info(
        "%s",
        tabulate(
            [[sg, f"{score:.4f}"] for sg, score in zip(space_groups, scores)],
            ["Space group", "score"],
        ),
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
