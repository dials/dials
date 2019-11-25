"""
A docstring

This can double as a helpful message which explains how the program is run.
"""

from __future__ import absolute_import, division, print_function

import logging
import sys

import libtbx.phil
import dials.util
import dials.util.log
from dials.util.options import OptionParser
from dials.util.options import flatten_reflections
from dials.util.options import flatten_experiments


# Define a logger
logger = logging.getLogger("dials.logger_name")

# Define the master PHIL scope for this program
phil_scope = libtbx.phil.parse(
    """
output {
    reflections = stronger.refl
        .type = path
    log = dials.command_name.log
        .type = path
}
bool_parameter = False
  .type = bool
integer_parameter = 0
  .type = int
"""
)


def do_stuff(experiments, reflections, params):
    """Write the behaviour of the program as functions and classes outside run()"""
    logger.info("Hello world!")
    # Include file output here
    logger.info("integer_parameter: %i", params.integer_parameter)
    logger.info("bool_parameter: %s", params.bool_parameter)
    reflections.as_file(params.output.reflections)


def run(args=None, phil=phil_scope):
    """
  run(args : [string] = sys.argv[1:],
      phil : libtbx.phil.scope = phil_scope)

  Run the script, parsing arguments found in 'args' (default: sys.argv[1:])
  and using the PHIL scope defined in 'phil' (default: phil_scope,
  the master PHIL scope for this program).

  Try to keep this function minimal, defining only what is necessary to run
  the program from the command line.
  """
    usage = "dials.command_name [options] experiments.json strong.pickle"

    parser = OptionParser(
        usage=usage,
        phil=phil,
        read_reflections=True,
        read_experiments=True,
        check_format=False,
        epilog=__doc__,
    )

    params, options = parser.parse_args(args=args, show_diff_phil=False)

    # Configure the logging
    dials.util.log.config(options.verbose, logfile=params.output.log)

    # Log the PHIL diff
    # Log the diff phil
    diff_phil = parser.diff_phil.as_str()
    if diff_phil:
        logger.info("The following parameters have been modified:\n%s", diff_phil)

    experiments = flatten_experiments(params.input.experiments)
    reflections = flatten_reflections(params.input.reflections)

    if len(reflections) != 1:
        sys.exit("Exactly one reflection file needed")
    if len(experiments) != 1:
        sys.exit("Exactly one 1 experiment list required")

    do_stuff(experiments, reflections[0], params)


if __name__ == "__main__":
    with dials.util.show_mail_on_error():
        run()
