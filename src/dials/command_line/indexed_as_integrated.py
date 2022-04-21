from __future__ import annotations

import logging

import iotbx.phil

import dials.util
from dials.array_family import flex
from dials.util import log
from dials.util.options import ArgumentParser, reflections_and_experiments_from_files
from dials.util.version import dials_version

help_message = """
Small modifications to an indexed spot list to allow it to be treated as if
it were integrated, for rapid feedback analysis. Do not use output data for
structure solution or refinement.
"""

phil_scope = iotbx.phil.parse(
    """
output {
  reflections = pseudo_integrated.refl
    .type = path
  log = "dials.indexed_as_integrated.log"
    .type = path
}
"""
)


logger = logging.getLogger("dials.command_line.indexed_as_integrated")


def indexed_as_integrated(reflections, params, experiments):
    """Small modifications to an indexed spot list to allow it to be
    treated as if it were integrated, for rapid feedback analysis."""

    # filter only indexed reflections & assign the summation integrated flag

    sel = reflections.get_flags(reflections.flags.indexed)
    reflections = reflections.select(sel)
    all = flex.bool(reflections.size(), True)
    reflections.set_flags(all, reflections.flags.integrated_sum)

    logger.info(f"\nSaved {sel.count(True)} of {sel.size()} reflections\n")

    # add resolution to reflections

    if "imageset_id" not in reflections:
        reflections["imageset_id"] = reflections["id"]

    reflections.centroid_px_to_mm(experiments)
    reflections.map_centroids_to_reciprocal_space(experiments)
    reflections["d"] = 1 / reflections["rlp"].norms()

    return reflections


@dials.util.show_mail_handle_errors()
def run(args=None):
    usage = "dials.indexed_as_integrated [options] indexed.refl indexed.expt"

    parser = ArgumentParser(
        usage=usage,
        phil=phil_scope,
        read_reflections=True,
        read_experiments=True,
        check_format=False,
        epilog=help_message,
    )

    params, _ = parser.parse_args(args, show_diff_phil=True)

    log.config(logfile=params.output.log)
    logger.info(
        dials_version()
        + "\n\n----------------------------------------------------------\n"
        "Do not use the output for structure solution or refinement"
        "\n----------------------------------------------------------"
    )

    reflections, experiments = reflections_and_experiments_from_files(
        params.input.reflections, params.input.experiments
    )

    if len(reflections) != 1 or len(experiments) != 1:
        parser.print_help()
        return

    reflections = reflections[0]

    pseudo_integrated = indexed_as_integrated(reflections, params, experiments)
    pseudo_integrated.as_file(params.output.reflections)

    logger.info(f"\nSaved pseudo-integrated data to: {params.output.reflections}\n")


if __name__ == "__main__":
    run()
