#!/usr/bin/env python
#
# dials.ssx_integrate.py
#
#  Copyright (C) 2021 Diamond Light Source
#
#  Author: James Beilsten-Edmands
#
#  This code is distributed under the BSD license, a copy of which is
#  included in the root directory of this package.


from __future__ import absolute_import, division

import concurrent.futures
import functools
import logging

import iotbx.phil
from dxtbx.model import ExperimentList
from libtbx import Auto
from libtbx.introspection import number_of_processors
from libtbx.utils import Sorry

from dials.algorithms.integration.ssx.potato_integrate import (
    PotatoIntegrator,
    PotatoOutputAggregator,
)
from dials.algorithms.integration.ssx.ssx_integrate import (
    OutputAggregator,
    generate_html_report,
)
from dials.algorithms.integration.ssx.stills_integrate import StillsIntegrator
from dials.array_family import flex
from dials.util import log, show_mail_handle_errors
from dials.util.options import OptionParser, flatten_experiments, flatten_reflections
from dials.util.version import dials_version

try:
    from typing import List
except ImportError:
    pass

logger = logging.getLogger("dials.ssx_integrate")

help_message = """

This program does profile modelling and integration for stills

"""

# Create the phil scope
phil_scope = iotbx.phil.parse(
    """
  algorithm = *potato stills
    .type = choice
  nproc=Auto
    .type = int
  output {
    batch_size = 50
      .type = int
      .help = "Number of images to save in each output file"
    log = "dials.ssx_integrate.log"
      .type = str
    html = "dials.ssx_integrate.html"
      .type = str
  }

  potato {
    include scope dials.algorithms.profile_model.potato.potato.phil_scope
  }

  stills {
    include scope dials.algorithms.integration.integrator.phil_scope
    include scope dials.algorithms.profile_model.factory.phil_scope
    include scope dials.algorithms.spot_prediction.reflection_predictor.phil_scope
    include scope dials.algorithms.integration.stills_significance_filter.phil_scope
    include scope dials.algorithms.integration.kapton_correction.absorption_phil_scope
  }

  debug {
    output {
      shoeboxes = False
        .type = bool
    }
  }

""",
    process_includes=True,
)

phil_overrides = phil_scope.fetch(
    source=iotbx.phil.parse(
        """\
stills {
    profile {
        gaussian_rs {
            min_spots {
                overall=0
            }
        }
        fitting = False
    }
    integration {
        background {
            simple {
                outlier {
                    algorithm = Null
                }
            }
        }
    }
}
potato {
    refinement {
        n_cycles = 1
    }
}
"""
    )
)
working_phil = phil_scope.fetch(sources=[phil_overrides])

loggers_to_disable = ["dials", "dials.array_family.flex_ext"]

loggers_to_disable_for_stills = loggers_to_disable + [
    "dials.algorithms.integration.integrator",
    "dials.algorithms.profile_model.gaussian_rs.calculator",
    "dials.command_line.integrate",
    "dials.algorithms.spot_prediction.reflection_predictor",
]


def disable_loggers(lognames: List[str]) -> None:
    for logname in lognames:
        logging.getLogger(logname).disabled = True


def process_one_image_potato_integrator(experiment, table, params):

    disable_loggers(loggers_to_disable)  # disable the loggers within each process

    integrator = PotatoIntegrator(params.potato, collect_data=params.output.html)
    try:
        experiment, table, collector = integrator.run(experiment, table)
    except RuntimeError as e:
        logger.info(f"Processing failed due to error: {e}")
        return (None, None, None)
    else:
        return experiment, table, collector


def process_one_image_stills_integrator(experiment, table, params):

    disable_loggers(
        loggers_to_disable_for_stills
    )  # disable the loggers within each process

    integrator = StillsIntegrator(params.stills, collect_data=params.output.html)
    try:
        experiment, table, collector = integrator.run(experiment, table)
    except RuntimeError as e:
        logger.info(f"Processing failed due to error: {e}")
        return (None, None, None)
    else:
        return experiment, table, collector


def setup(reflections, params):
    # calculate the batches for processing
    batches = list(range(0, len(reflections), params.output.batch_size))
    batches.append(len(reflections))

    # Note, memory processing logic can go here
    if params.nproc is Auto:
        params.nproc = number_of_processors(return_value_if_unknown=1)
    logger.info(f"Using {params.nproc} processes for integration")

    # aggregate some output for json, html etc
    if params.algorithm == "potato":
        process = process_one_image_potato_integrator
        aggregator = PotatoOutputAggregator()
    elif params.algorithm == "stills":
        process = process_one_image_stills_integrator
        aggregator = OutputAggregator()
    else:
        raise ValueError("Invalid algorithm choice")

    configuration = {
        "process": process,
        "aggregator": aggregator,
        "params": params,
    }

    return batches, configuration


def process_batch(sub_tables, sub_expts, configuration, batch_offset=0):
    integrated_reflections = flex.reflection_table()
    with concurrent.futures.ProcessPoolExecutor(
        max_workers=configuration["params"].nproc
    ) as pool:
        futures = {
            pool.submit(
                configuration["process"], expt, table, configuration["params"]
            ): i
            for i, (table, expt) in enumerate(zip(sub_tables, sub_expts))
        }
        tables_list = [0] * len(sub_expts)
        expts_list = [0] * len(sub_expts)
        for future in concurrent.futures.as_completed(futures):
            try:
                expt, refls, collector = future.result()
                j = futures[future]
            except Exception as e:
                logger.info(e)
            else:
                if refls and expt:
                    logger.info(f"Processed image {j+batch_offset+1}")
                    tables_list[j] = refls
                    expts_list[j] = expt
                    configuration["aggregator"].add_dataset(
                        collector, j + batch_offset + 1
                    )

        expts_list = list(filter(lambda a: a != 0, expts_list))
        integrated_experiments = ExperimentList(expts_list)

        n_integrated = 0
        for _ in range(len(tables_list)):
            table = tables_list.pop(0)
            if not table:
                continue
            # renumber actual id before extending
            ids_map = dict(table.experiment_identifiers())
            assert len(ids_map) == 1, ids_map
            del table.experiment_identifiers()[list(ids_map.keys())[0]]
            table["id"] = flex.int(table.size(), n_integrated)
            table.experiment_identifiers()[n_integrated] = list(ids_map.values())[0]
            n_integrated += 1
            if not configuration["params"].debug.output.shoeboxes:
                del table["shoebox"]
            integrated_reflections.extend(table)
            del table

        integrated_reflections.assert_experiment_identifiers_are_consistent(
            integrated_experiments
        )
    return integrated_experiments, integrated_reflections


@show_mail_handle_errors()
def run(args: List[str] = None, phil=working_phil) -> None:
    """Run from the command-line."""
    usage = ""

    parser = OptionParser(
        usage=usage,
        phil=phil,
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

    batches, configuration = setup(reflections, params)

    # determine suitable output filenames
    template = "{prefix}_{index:0{maxindexlength:d}d}.{extension}"
    experiments_template = functools.partial(
        template.format,
        prefix="integrated",
        maxindexlength=len(str(len(batches) - 1)),
        extension="expt",
    )
    reflections_template = functools.partial(
        template.format,
        prefix="integrated",
        maxindexlength=len(str(len(batches) - 1)),
        extension="refl",
    )

    # now process each batch, and do parallel processing within a batch
    for i, b in enumerate(batches[:-1]):
        end_ = batches[i + 1]
        logger.info(f"Processing images {b+1} to {end_}")
        sub_tables = reflections[b:end_]
        sub_expts = experiments[b:end_]

        integrated_experiments, integrated_reflections = process_batch(
            sub_tables, sub_expts, configuration, batch_offset=b
        )

        experiments_filename = experiments_template(index=i)
        reflections_filename = reflections_template(index=i)
        # Save the reflections
        logger.info(
            f"Saving {integrated_reflections.size()} reflections to {reflections_filename}"
        )
        integrated_reflections.as_file(reflections_filename)
        logger.info(f"Saving the experiments to {experiments_filename}")
        integrated_experiments.as_file(experiments_filename)

    # now generate a html report using the aggregated data.
    plots = configuration["aggregator"].make_plots()

    if params.output.html and plots:
        logger.info(f"Writing html report to {params.output.html}")
        generate_html_report(plots, params.output.html)


if __name__ == "__main__":
    run()
