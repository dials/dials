#!/usr/bin/env python
#
# dials.potato.py
#
#  Copyright (C) 2018 Diamond Light Source
#
#  Author: James Parkhurst
#
#  This code is distributed under the BSD license, a copy of which is
#  included in the root directory of this package.

# LIBTBX_SET_DISPATCHER_NAME dials.potato

from __future__ import absolute_import, division

import logging

from libtbx import phil
from libtbx.phil import parse
from libtbx.utils import Sorry

from dials.algorithms.profile_model.potato.potato import Integrator
from dials.array_family import flex
from dials.util import log, show_mail_handle_errors
from dials.util.options import OptionParser, flatten_experiments, flatten_reflections
from dials.util.version import dials_version

try:
    from typing import List
except ImportError:
    pass

logger = logging.getLogger("dials")

help_message = """

This script does profile modelling for stills

"""

# Create the phil scope
phil_scope = parse(
    """

  output {
    experiments = "integrated.expt"
      .type = str
      .help = "The output experiments"

    reflections = "integrated.refl"
      .type = str
      .help = "The output reflections"

    log = "dials.potato.log"
      .type = str
  }

  include scope dials.algorithms.profile_model.potato.potato.phil_scope

""",
    process_includes=True,
)


@show_mail_handle_errors()
def run(args: List[str] = None, phil: phil.scope = phil_scope) -> None:
    """Run from the command-line."""
    usage = ""

    parser = OptionParser(
        usage=usage,
        phil=phil_scope,
        epilog=help_message,
        read_experiments=True,
        read_reflections=True,
    )
    # Check the number of arguments is correct

    # Parse the command line
    params, options = parser.parse_args(show_diff_phil=False)
    reflections = flatten_reflections(params.input.reflections)
    experiments = flatten_experiments(params.input.experiments)
    if len(reflections) == 0 or len(experiments) == 0:
        parser.print_help()
        return
    elif len(reflections) != 1:
        raise Sorry("more than 1 reflection file was given")
    elif len(experiments) == 0:
        raise Sorry("no experiment list was specified")
    reflections = reflections[0]

    # Configure logging
    log.config(verbosity=options.verbose, logfile=params.output.log)
    logger.info(dials_version())

    # Log the diff phil
    diff_phil = parser.diff_phil.as_str()
    if diff_phil != "":
        logger.info("The following parameters have been modified:\n")
        logger.info(diff_phil)

    # Contruct the integrator
    reflections["id"] = flex.int(reflections.size(), 0)
    integrator = Integrator(experiments, reflections, params)

    # Do cycles of indexing and refinement
    for i in range(params.refinement.n_macro_cycles):
        integrator.reindex_strong_spots()
        integrator.integrate_strong_spots()
        integrator.refine()

    # Do the integration
    integrator.predict()
    integrator.integrate()

    # Get the reflections
    reflections = integrator.reflections
    experiments = integrator.experiments

    # Save the reflections
    logger.info(
        f"Saving {reflections.size()} reflections to {params.output.reflections}"
    )
    reflections.as_file(params.output.reflections)
    logger.info(f"Saving the experiments to {params.output.experiments}")
    experiments.as_file(params.output.experiments)


if __name__ == "__main__":
    run()
