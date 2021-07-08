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

from dxtbx.model import ExperimentList
from libtbx import phil
from libtbx.phil import parse
from libtbx.utils import Sorry

from dials.algorithms.profile_model.potato.potato import (
    Integrator,
    generate_html_report,
)
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
    html = "dials.potato.html"
      .type = str
  }

  include scope dials.algorithms.profile_model.potato.potato.phil_scope

""",
    process_includes=True,
)


def process_one_still(experiment, table, params):
    single_elist = ExperimentList([experiment])
    ids_map = dict(table.experiment_identifiers())
    table["id"] = flex.int(table.size(), 0)  # set to zero so can integrate (
    # this is how integration finds the image in the imageset)
    del table.experiment_identifiers()[list(ids_map.keys())[0]]
    table.experiment_identifiers()[0] = list(ids_map.values())[0]

    integrator = Integrator(single_elist, table, params)
    # Do cycles of indexing and refinement
    for i in range(params.refinement.n_macro_cycles):
        try:
            integrator.reindex_strong_spots()
        except RuntimeError as e:
            logger.info(f"Processing failed due to error: {e}")
            return (None, None, None)
        else:
            integrator.integrate_strong_spots()
            try:
                integrator.refine()
            except RuntimeError as e:
                logger.info(f"Processing failed due to error: {e}")
                return (None, None, None)

    # Do the integration
    integrator.predict()
    integrator.integrate()

    # Get the reflections
    return integrator.experiments, integrator.reflections, integrator.plots_data


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

    if len(reflections) != 1:
        raise Sorry(
            "Only a single reflection table file can be input (this can be a multi-still table)"
        )

    # Configure logging
    log.config(verbosity=options.verbose, logfile=params.output.log)
    logger.info(dials_version())

    # Log the diff phil
    diff_phil = parser.diff_phil.as_str()
    if diff_phil != "":
        logger.info("The following parameters have been modified:\n")
        logger.info(diff_phil)

    ## FIXME - experiment identifiers approach wont work if input strong.refl and refined.expt
    # for now - check image path and update identifiers to that of refined.expt?
    if len(set(reflections[0]["id"]).difference({-1})) > 1:
        logger.info("Attempting to split multi-still reflection table")
        reflections = reflections[0].split_by_experiment_id()
        if not (len(reflections) == len(experiments)):
            raise Sorry(
                "Unequal number of reflection tables and experiments after splitting"
            )

    integrated_experiments = ExperimentList()
    integrated_reflections = flex.reflection_table()

    n_integrated = 0
    for n, (table, expt) in enumerate(zip(reflections, experiments)):
        experiments, reflections, plots_data = process_one_still(expt, table, params)
        if not (experiments and reflections):
            continue
        # renumber actual id before extending
        ids_map = dict(reflections.experiment_identifiers())
        assert len(ids_map) == 1
        del reflections.experiment_identifiers()[list(ids_map.keys())[0]]
        reflections["id"] = flex.int(reflections.size(), n_integrated)
        reflections.experiment_identifiers()[n_integrated] = list(ids_map.values())[0]
        n_integrated += 1
        integrated_experiments.extend(experiments)
        integrated_reflections.extend(reflections)

    integrated_reflections.assert_experiment_identifiers_are_consistent(
        integrated_experiments
    )

    # Save the reflections
    logger.info(
        f"Saving {integrated_reflections.size()} reflections to {params.output.reflections}"
    )
    integrated_reflections.as_file(params.output.reflections)
    logger.info(f"Saving the experiments to {params.output.experiments}")
    integrated_experiments.as_file(params.output.experiments)

    if params.output.html and plots_data:
        # FIXME currently only plots for last image processed
        generate_html_report(plots_data, params.output.html)


if __name__ == "__main__":
    run()
