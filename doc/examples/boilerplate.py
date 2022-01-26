"""
A docstring

This can double as a helpful message which explains how the program is run.
"""

from __future__ import annotations

import logging
import sys

# We need to parse command-line arguments to PHIL scopes.
import libtbx.phil

# Information about the experiment geometry and meta-data are recorded in
# dxtbx.model.Experiment objects, collated in ExperimentList objects.
from dxtbx.model import ExperimentList

# All command-line DIALS programs should run with dials.util.show_mail_handle_errors.
import dials.util

# The logging module is used to raise log messages.  Additionally, dials.util.log
# sets up handlers and filters for a consistent logging style across DIALS.
import dials.util.log

# Often, we deal with flex arrays and reflection tables to handle reflection data.
from dials.array_family import flex

# The DIALS option parser is based on the (old) standard Python option parser,
# but contains customisations such as the parsing of PHIL parameters.
# flatten_experiments & flatten_reflections are useful for combining multiple input
# experiment lists and reflection tables into a single instance of each.
from dials.util.options import ArgumentParser, flatten_experiments, flatten_reflections

# Useful to know what version of DIALS we are running
from dials.util.version import dials_version

try:
    from typing import List
except ImportError:
    pass


# Define a logger.
logger = logging.getLogger("dials.boilerplate")

# Define the master PHIL scope for this program.
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


def do_boilerplate(
    experiments: ExperimentList,
    reflections: flex.reflection_table,
    params: libtbx.phil.scope_extract,
) -> (ExperimentList, flex.reflection_table):
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
        RuntimeError:  if someone says 'goose'.
    """
    logger.info("Hello world!")

    # Here's an example of an error that might be raised, as documented above.
    if "goose" in ["duck", "duck", "duck"]:
        raise RuntimeError("Quick, run!")

    logger.info("The input reflection table contains %d reflections.", len(reflections))
    logger.info(
        "The input experiment list contains %d imagesets.", len(experiments.imagesets())
    )

    logger.info("integer_parameter: %i", params.integer_parameter)
    logger.info("bool_parameter: %s", params.bool_parameter)

    return experiments, reflections


@dials.util.show_mail_on_error()
def run(args: List[str] = None, phil: libtbx.phil.scope = phil_scope) -> None:
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

    parser = ArgumentParser(
        usage=usage,
        phil=phil,
        read_reflections=True,
        read_experiments=True,
        check_format=False,
        epilog=__doc__,
    )

    params, options = parser.parse_args(args=args, show_diff_phil=False)

    # Configure the logging.
    dials.util.log.config(options.verbose, logfile=params.output.log)

    # Log the dials version
    logger.info(dials_version())

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
        sys.exit("Exactly one reflection file needed.")
    if len(experiments) != 1:
        sys.exit("Exactly one experiment list required.")

    # Do whatever this program is supposed to do.
    experiments, reflections = do_boilerplate(experiments, reflections[0], params)

    # Do the file output here.
    logger.info("Writing the reflection table to %s", params.output.reflections)
    reflections.as_file(params.output.reflections)


# Keep this minimal. Calling run() should do exactly the same thing as running this
if __name__ == "__main__":
    run()
