from __future__ import annotations

import logging
import random
import sys

import numpy as np
from jinja2 import ChoiceLoader, Environment, PackageLoader

import iotbx.phil

from dials.algorithms.correlation.analysis import CorrelationMatrix
from dials.array_family import flex
from dials.util import log, show_mail_handle_errors
from dials.util.multi_dataset_handling import (
    assign_unique_identifiers,
    parse_multiple_datasets,
)
from dials.util.options import ArgumentParser, reflections_and_experiments_from_files
from dials.util.version import dials_version

logger = logging.getLogger("dials.algorithms.correlation.analysis")

phil_scope = iotbx.phil.parse(
    """\
include scope dials.algorithms.correlation.analysis.working_phil

output {
  log = dials.correlation_matrix.log
    .type = str
    .help = "The log name"
  html = dials.correlation_matrix.html
    .type = path
    .help = "Filename for the html report"
  json = None
    .type = str
    .help = "Filename for the cluster information output in json format"
}
significant_clusters {
  output = False
    .type = bool
    .help = "Toggle to output expt/refl files for significant clusters as determined by OPTICS clustering on cosine angle coordinates"
}
""",
    process_includes=True,
)


help_message = """
This module implements a subset of methods used in dials.cosym to perform
correlation and cosine similarity based clustering methods. Data should be passed
through dials.cosym first to implement consistent symmetry.

Examples::

  dials.correlation_matrix symmetrized.expt symmetrized.refl

"""


@show_mail_handle_errors()
def run(args=None):
    usage = "dials.correlation_matrix [options] symmetrized.expt symmetrized.refl"

    parser = ArgumentParser(
        usage=usage,
        phil=phil_scope,
        read_reflections=True,
        read_experiments=True,
        check_format=False,
        epilog=help_message,
    )

    params, options, args = parser.parse_args(
        args=args, show_diff_phil=False, return_unhandled=True
    )

    # Configure the logging
    log.config(verbosity=options.verbose, logfile=params.output.log)

    logger.info(dials_version())

    # Log the diff phil
    diff_phil = parser.diff_phil.as_str()
    if diff_phil != "":
        logger.info("The following parameters have been modified:\n")
        logger.info(diff_phil)

    if params.seed is not None:
        flex.set_random_seed(params.seed)
        np.random.seed(params.seed)
        random.seed(params.seed)

    if not params.input.experiments or not params.input.reflections:
        parser.print_help()
        sys.exit()

    reflections, experiments = reflections_and_experiments_from_files(
        params.input.reflections, params.input.experiments
    )

    reflections = parse_multiple_datasets(reflections)
    if len(experiments) != len(reflections):
        sys.exit(
            "Mismatched number of experiments and reflection tables found: %s & %s."
            % (len(experiments), len(reflections))
        )
    if len(experiments) < 2:
        sys.exit(
            "At least 2 datasets are needed for cluster analysis. Please re-run with more datasets."
        )
    try:
        experiments, reflections = assign_unique_identifiers(experiments, reflections)
        matrices = CorrelationMatrix(
            experiments=experiments, reflections=reflections, params=params
        )
    except ValueError as e:
        sys.exit(e)

    matrices.calculate_matrices()

    if params.significant_clusters.output:
        matrices.output_clusters()
    else:
        logger.info(
            "For separated clusters in DIALS .expt/.refl output please re-run with significant_clusters.output=True"
        )

    if params.output.json:
        matrices.output_json()

    if params.output.html:
        matrices.convert_to_html_json()

        loader = ChoiceLoader(
            [PackageLoader("dials", "templates"), PackageLoader("dials", "templates")]
        )
        env = Environment(loader=loader)

        template = env.get_template("clusters.html")
        html = template.stream(
            page_title="DIALS Correlation Matrix",
            cc_cluster_json=matrices.cc_json,
            cc_cluster_table=matrices.cc_table,
            cos_angle_cluster_json=matrices.cos_json,
            cos_angle_cluster_table=matrices.cos_table,
            image_range_tables=[matrices.table_list],
            cosym_graphs=matrices.rij_graphs,
        )

        logger.info(
            f"Saving graphical output of correlation matrices to {params.output.html}."
        )
        html.dump(params.output.html, errors="xmlcharrefreplace")


if __name__ == "__main__":
    run()
