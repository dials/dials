"""
A docstring

This can double as a helpful message which explains how the program is run.
"""

from __future__ import absolute_import, division, print_function

import logging
import sys
from typing import List  # noqa: F401  Flake8 doesn't recognise Python 2 typing.

import libtbx.phil
from dials.array_family import flex
import dials.util
import dials.util.log
from dials.util.options import OptionParser, flatten_reflections, flatten_experiments
from dxtbx.model import ExperimentList


# Define a logger
logger = logging.getLogger("dials.logger_name")

# Define the master PHIL scope for this program
phil_scope = libtbx.phil.parse(
    """
    output {
        reflections = stronger.refl
            .type = path
            .help = "Help strings are helpful, wherever possible."
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
    # type: (ExperimentList, flex.reflection_table, libtbx.phil.scope_extract) -> None
    """
    Write the behaviour of the program as functions and classes outside run().

    Don't include file output here, remember that this function may be re-used
    elsewhere by someone who doesn't need the output written immediately to file.

    It can be especially helpful to document any expected exceptions that might be
    raised, in order to keep track of what needs to be handled in any code that
    re-uses this function.

    Args:
        experiments:  An experiment list.
        reflections:  A reflection table.
        params:       Some parameters, in the form of a scope_extract object,
                      which is the usable form of a parsed PHIL scope.

    Raises:
        RuntimeError:  if you pass an empty experiment list.
    """
    logger.info("Hello world!")

    # Here's an example of an error that might be raised, as documented above.
    if not experiments:
        raise RuntimeError

    logger.info("The input reflection table contains %d reflections.", len(reflections))
    logger.info(
        "The input experiment list contains %d imagesets.", len(experiments.imagesets)
    )

    logger.info("integer_parameter: %i", params.integer_parameter)
    logger.info("bool_parameter: %s", params.bool_parameter)


def run(args=None, phil=phil_scope):  # type: (List[str], libtbx.phil.scope) -> None
    """
    Check command-line input and call other functions to do the legwork.

    Run the script, parsing arguments found in 'args' and using the PHIL scope
    defined in 'phil'.

    Try to keep this function minimal, defining only what is necessary to run
    the program from the command line.

    Args:
        args: The arguments supplied by the user (default: sys.argv[1:])
        phil: The PHIL scope definition (default: phil_scope, the master PHIL scope
        for this program).
    """
    usage = "dials.command_name [options] imported.expt strong.refl"

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

    # Log the difference between the PHIL scope definition and the active PHIL scope,
    # which will include the parsed user inputs.
    diff_phil = parser.diff_phil.as_str()
    if diff_phil:
        logger.info("The following parameters have been modified:\n%s", diff_phil)

    # These functions are commonly used to collate the input.
    experiments = flatten_experiments(params.input.experiments)
    reflections = flatten_reflections(params.input.reflections)

    # You might well wish to check here that the command-line input is appropriate.
    if len(reflections) != 1:
        sys.exit("Exactly one reflection file needed")
    if len(experiments) != 1:
        sys.exit("Exactly one 1 experiment list required")

    # Do whatever this program is supposed to do.
    do_stuff(experiments, reflections[0], params)

    # Do the file output here.
    reflections.as_file(params.output.reflections)


# Keep this minimal.  Try to keep the command-line behaviour neatly encapsulated in run.
if __name__ == "__main__":
    with dials.util.show_mail_on_error():
        run()
