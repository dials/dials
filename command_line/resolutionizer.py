from __future__ import absolute_import, division, print_function

import logging
import sys

import libtbx.phil

from dials.util import resolutionizer
from dials.util import log
from dials.util.options import OptionParser, reflections_and_experiments_from_files
from dials.util.version import dials_version
from dials.util.multi_dataset_handling import parse_multiple_datasets

logger = logging.getLogger("dials.resolutionizer")

help_message = """
"""


phil_scope = libtbx.phil.parse(
    """
include scope dials.util.resolutionizer.phil_defaults

output {
  log = dials.resolutionizer.log
    .type = path
}
""",
    process_includes=True,
)


def run(args):
    usage = (
        "dials.resolutionizer [options] (scaled.expt scaled.refl | scaled_unmerged.mtz)"
    )

    parser = OptionParser(
        usage=usage,
        phil=phil_scope,
        read_reflections=True,
        read_experiments=True,
        check_format=False,
        epilog=help_message,
    )

    params, options, unhandled = parser.parse_args(
        return_unhandled=True, show_diff_phil=True
    )

    reflections, experiments = reflections_and_experiments_from_files(
        params.input.reflections, params.input.experiments
    )
    if (not reflections or not experiments) and not unhandled:
        parser.print_help()
        return

    if reflections and experiments and unhandled:
        sys.exit(
            "Must provide either scaled unmerged mtz OR dials-format scaled reflections and experiments files"
        )

    # Configure the logging
    log.config(logfile=params.output.log, verbosity=options.verbose)
    logger.info(dials_version())

    if len(unhandled) == 1:
        scaled_unmerged = unhandled[0]
        m = resolutionizer.Resolutionizer.from_unmerged_mtz(
            scaled_unmerged, params.resolutionizer
        )
    else:
        reflections = parse_multiple_datasets(reflections)
        m = resolutionizer.Resolutionizer.from_reflections_and_experiments(
            reflections, experiments, params.resolutionizer
        )

    m.resolution_auto()


if __name__ == "__main__":
    run(sys.argv[1:])
