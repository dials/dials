# LIBTBX_SET_DISPATCHER_NAME dev.dials.ssx_integrate
#!/usr/bin/env python
#
# dials.ssx_integrate.py
#
#  Copyright (C) 2022 Diamond Light Source
#
#  Author: James Beilsten-Edmands
#
#  This code is distributed under the BSD license, a copy of which is
#  included in the root directory of this package.
"""
This program rums profile modelling and integration on indexed results from a
still sequence i.e. SSX data. This scripts uses parts of the regular DIALS
integration code, using either the ellipsoid or stills integrator algorithms.

The ellipsoid algorithm refines the unit cell, orientation and a 3D ellipsoidal
mosaicity parameterisation for each crystal, by assessing the pixel-intensity
distribution of the strong spots. The integrated data are saved in batches to
hep with memory management. A html output report is generated, showing integration
and clutering statistics.

Further program documentation can be found at dials.github.io/ssx_processing_guide.html

Usage:
    dev.dials.ssx_integrate indexed.expt indexed.refl
    dev.dials.ssx_integrate refined.expt refined.refl
    dev.dials.ssx_integrate indexed.expt indexed.refl algorithm=stills
"""

from __future__ import annotations

import concurrent.futures
import copy
import json
import logging

import iotbx.phil
from cctbx import crystal
from dxtbx.model import ExperimentList
from libtbx import Auto
from libtbx.introspection import number_of_processors
from libtbx.utils import Sorry

from dials.algorithms.indexing.ssx.analysis import report_on_crystal_clusters
from dials.algorithms.integration.ssx.ellipsoid_integrate import (
    EllipsoidIntegrator,
    EllipsoidOutputAggregator,
)
from dials.algorithms.integration.ssx.ssx_integrate import (
    OutputAggregator,
    generate_html_report,
)
from dials.algorithms.integration.ssx.stills_integrate import StillsIntegrator
from dials.array_family import flex
from dials.util import log, show_mail_handle_errors
from dials.util.options import ArgumentParser, flatten_experiments, flatten_reflections
from dials.util.version import dials_version

try:
    from typing import List
except ImportError:
    pass

logger = logging.getLogger("dials.ssx_integrate")

# Create the phil scope
phil_scope = iotbx.phil.parse(
    """
  algorithm = *ellipsoid stills
    .type = choice
  nproc=Auto
    .type = int
  image_range = None
    .type = str
    .help = "Only process this image range. Input in the format start:end."
            "This range is an inclusive range."
  output {
    batch_size = 50
      .type = int
      .help = "Number of images to save in each output file"
    log = "dials.ssx_integrate.log"
      .type = str
    html = "dials.ssx_integrate.html"
      .type = str
    json = None
      .type = str
  }

  ellipsoid {
    include scope dials.algorithms.profile_model.ellipsoid.algorithm.ellipsoid_algorithm_phil_scope
  }

  include scope dials.algorithms.integration.integrator.phil_scope
  include scope dials.algorithms.profile_model.factory.phil_scope
  include scope dials.algorithms.spot_prediction.reflection_predictor.phil_scope
  include scope dials.algorithms.integration.stills_significance_filter.phil_scope
  include scope dials.algorithms.integration.kapton_correction.absorption_phil_scope


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
profile {
    gaussian_rs {
        min_spots {
            overall=0
        }
    }
    fitting = False
    ellipsoid {
        refinement {
            n_cycles = 1
        }
    }
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
"""
    )
)
working_phil = phil_scope.fetch(sources=[phil_overrides])

working_phil.adopt_scope(
    iotbx.phil.parse(
        """
    individual_log_verbosity = 1
    .type =int
"""
    )
)

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


def process_one_image_ellipsoid_integrator(experiment, table, params):

    if params.individual_log_verbosity < 2:
        disable_loggers(loggers_to_disable)  # disable the loggers within each process
    elif params.individual_log_verbosity == 2:
        for name in loggers_to_disable:
            logging.getLogger(name).setLevel(logging.INFO)

    collect_data = params.output.html or params.output.json
    integrator = EllipsoidIntegrator(params, collect_data)
    try:
        experiment, table, collector = integrator.run(experiment, table)
    except RuntimeError as e:
        logger.info(f"Processing failed due to error: {e}")
        return (None, None, None)
    else:
        return experiment, table, collector


def process_one_image_stills_integrator(experiment, table, params):

    if params.individual_log_verbosity < 2:
        disable_loggers(
            loggers_to_disable_for_stills
        )  # disable the loggers within each process
    elif params.individual_log_verbosity == 2:
        for name in loggers_to_disable_for_stills:
            logging.getLogger(name).setLevel(logging.INFO)

    collect_data = params.output.html or params.output.json
    integrator = StillsIntegrator(params, collect_data)
    try:
        experiment, table, collector = integrator.run(experiment, table)
    except RuntimeError as e:
        logger.info(f"Processing failed due to error: {e}")
        return (None, None, None)
    else:
        return experiment, table, collector


def setup(reflections, params):
    # calculate the batches for processing
    if params.image_range:
        if not len(params.image_range.split(":")) == 2:
            raise ValueError("Image range must be given in the form first:last")
        first, last = params.image_range.split(":")
        try:
            first = int(first)
        except ValueError:
            raise ValueError(f"Issue interpreting {first} as an integer")
        try:
            last = int(last)
        except ValueError:
            raise ValueError(f"Issue interpreting {last} as an integer")
        lowest_index = max(0, first - 1)
        highest_image = min(last, len(reflections))
        batches = list(range(lowest_index, highest_image, params.output.batch_size))
        batches.append(highest_image)
    else:
        batches = list(range(0, len(reflections), params.output.batch_size))
        batches.append(len(reflections))

    # Note, memory processing logic can go here
    if params.nproc is Auto:
        params.nproc = number_of_processors(return_value_if_unknown=1)
    logger.info(f"Using {params.nproc} processes for integration")

    # aggregate some output for json, html etc
    if params.algorithm == "ellipsoid":
        process = process_one_image_ellipsoid_integrator
        aggregator = EllipsoidOutputAggregator()
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


def run_integration(reflections, experiments, params):
    assert len(reflections) == len(experiments)
    batches, configuration = setup(reflections, params)

    # now process each batch, and do parallel processing within a batch
    for i, b in enumerate(batches[:-1]):
        end_ = batches[i + 1]
        logger.info(f"Processing images {b+1} to {end_}")
        sub_tables = reflections[b:end_]
        sub_expts = experiments[b:end_]

        integrated_experiments, integrated_reflections = process_batch(
            sub_tables, sub_expts, configuration, batch_offset=b
        )
        yield integrated_experiments, integrated_reflections, configuration[
            "aggregator"
        ]


@show_mail_handle_errors()
def run(args: List[str] = None, phil=working_phil) -> None:
    """
    Run dev.dials.ssx_integrate from the command-line.

    This program takes an indexed experiment list and reflection table and
    performs parallelised integration for synchrotron serial crystallography
    experiments. The programs acts as a wrapper to run one of two algorithms,
    the stills integrator or the 'ellipsoid' integrator (which uses a generalised
    ellipsoidal profile model). Analysis statistics are captured and output as
    a html report, while the output data are saved in batches for memory
    management.
    """

    parser = ArgumentParser(
        usage="dev.dials.ssx_integrate indexed.expt indexed.refl [options]",
        phil=phil,
        epilog=__doc__,
        read_experiments=True,
        read_reflections=True,
    )
    # Check the number of arguments is correct

    # Parse the command line
    params, options = parser.parse_args(args=args, show_diff_phil=False)
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
    params.individual_log_verbosity = options.verbose
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

    integrated_crystal_symmetries = []

    for i, (int_expt, int_refl, aggregator) in enumerate(
        run_integration(reflections, experiments, params)
    ):

        reflections_filename = f"integrated_{i+1}.refl"
        experiments_filename = f"integrated_{i+1}.expt"
        logger.info(f"Saving {int_refl.size()} reflections to {reflections_filename}")
        int_refl.as_file(reflections_filename)
        logger.info(f"Saving the experiments to {experiments_filename}")
        int_expt.as_file(experiments_filename)

        integrated_crystal_symmetries.extend(
            [
                crystal.symmetry(
                    unit_cell=copy.deepcopy(cryst.get_unit_cell()),
                    space_group=copy.deepcopy(cryst.get_space_group()),
                )
                for cryst in int_expt.crystals()
            ]
        )

    plots = {}
    cluster_plots, _ = report_on_crystal_clusters(
        integrated_crystal_symmetries,
        make_plots=(params.output.html or params.output.json),
    )

    if params.output.html or params.output.json:
        # now generate plots using the aggregated data.
        plots = aggregator.make_plots()
        plots.update(cluster_plots)

    if params.output.html and plots:
        logger.info(f"Writing html report to {params.output.html}")
        generate_html_report(plots, params.output.html)
    if params.output.json and plots:
        logger.info(f"Saving plot data in json format to {params.output.json}")
        with open(params.output.json, "w") as outfile:
            json.dump(plots, outfile, indent=2)

    logger.info(
        "Further program documentation can be found at dials.github.io/ssx_processing_guide.html"
    )


if __name__ == "__main__":
    run()
