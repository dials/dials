# LIBTBX_SET_DISPATCHER_NAME dials.resolutionizer
# LIBTBX_SET_DISPATCHER_NAME xia2.resolutionizer

import logging
import sys

from dials.util import resolutionizer
from dials.util import log
from dials.util.version import dials_version


logger = logging.getLogger("dials.resolutionizer")


def run(args):
    working_phil = resolutionizer.phil_defaults
    interp = working_phil.command_line_argument_interpreter(home_scope="resolutionizer")
    params, unhandled = interp.process_and_fetch(
        args, custom_processor="collect_remaining"
    )
    params = params.extract().resolutionizer
    if len(unhandled) == 0:
        working_phil.show()
        exit()

    # Configure the logging
    log.config(logfile="dials.resolutionizer.log")
    logger.info(dials_version())

    assert len(unhandled) == 1
    scaled_unmerged = unhandled[0]

    m = resolutionizer.Resolutionizer.from_unmerged_mtz(scaled_unmerged, params)
    m.resolution_auto()


if __name__ == "__main__":
    run(sys.argv[1:])
