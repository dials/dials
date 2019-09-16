# LIBTBX_SET_DISPATCHER_NAME dials.resolutionizer
# LIBTBX_SET_DISPATCHER_NAME xia2.resolutionizer

import logging
import sys

from dials.util import resolutionizer
from dials.util import log
from dials.util.options import OptionParser
from dials.util.options import flatten_reflections, flatten_experiments
from dials.util.version import dials_version

logger = logging.getLogger("dials.resolutionizer")

help_message = """
"""


phil_scope = """
include scope dials.util.resolutionizer.phil_defaults
"""


def run(args):
    usage = (
        "dials.resolutionizer [options] scaled.expt scaled.refl | scaled_unmerged.mtz"
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

    # Configure the logging
    log.config(logfile=params.output.log)
    logger.info(dials_version())

    reflections = flatten_reflections(params.input.reflections)
    experiments = flatten_experiments(params.input.experiments)
    if len(reflections) == 0 or len(experiments) == 0 or len(unhandled) == 1:
        parser.print_help()
        return

    # Configure the logging
    log.config(logfile="dials.resolutionizer.log")
    logger.info(dials_version())

    if len(unhandled) == 1:
        scaled_unmerged = unhandled[0]
        m = resolutionizer.Resolutionizer.from_unmerged_mtz(scaled_unmerged, params)
    else:
        m = resolutionizer.Resolutionizer.from_reflections_and_experiments(
            reflections, experiments, params
        )

    m.resolution_auto()


if __name__ == "__main__":
    run(sys.argv[1:])
