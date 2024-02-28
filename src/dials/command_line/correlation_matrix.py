from __future__ import annotations

import logging
import random
import sys

import numpy as np
from jinja2 import ChoiceLoader, Environment, PackageLoader

import iotbx.phil

from dials.algorithms.correlation.analysis import CorrelationMatrix
from dials.array_family import flex
from dials.util import Sorry, log, show_mail_handle_errors
from dials.util.multi_dataset_handling import (
    assign_unique_identifiers,
    parse_multiple_datasets,
)
from dials.util.options import ArgumentParser, reflections_and_experiments_from_files
from dials.util.version import dials_version

logger = logging.getLogger("dials.algorithms.correlation.analysis")

phil_scope = iotbx.phil.parse(
    """\
include scope dials.algorithms.correlation.analysis.phil_scope

seed = 42
  .type = int(value_min=0)

output {
  log = dials.correlation_matrix.log
    .type = str
  html = dials.correlation_matrix.html
    .type = path
  cc_json = dials.correlation_matrix_cc.json
    .type = str
  cos_json = dials.correlation_matrix_cos.json
    .type = str
}
""",
    process_includes=True,
)


help_message = """
This module will implement a subset of methods used in dials.cosym to perform
correlation and cos angle based clustering methods. Data should be passed
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
        raise Sorry(
            "Mismatched number of experiments and reflection tables found: %s & %s."
            % (len(experiments), len(reflections))
        )
    if len(experiments) == 1:
        raise Sorry(
            "At least 2 datasets are needed for cluster analysis. Please re-run with more datasets."
        )
    try:
        experiments, reflections = assign_unique_identifiers(experiments, reflections)
        matrices = CorrelationMatrix(
            experiments=experiments, reflections=reflections, params=params
        )
    except ValueError as e:
        raise Sorry(e)

    matrices.calculate_matrices()
    matrices.output_json()

    loader = ChoiceLoader(
        [PackageLoader("dials", "templates"), PackageLoader("dials", "templates")]
    )
    env = Environment(loader=loader)

    template = env.get_template("clusters.html")
    html = template.render(
        page_title="DIALS Correlation Matrix",
        cc_cluster_json=matrices.cc_json,
        cos_angle_cluster_json=matrices.cos_json,
        image_range_tables=[matrices.table_list],
        cosym_graphs=matrices.rij_graphs,
    )

    logger.info(
        f"Saving graphical output of correlation matrices to {params.output.html}."
    )

    with open(params.output.html, "wb") as f:
        f.write(html.encode("utf-8", "xmlcharrefreplace"))


if __name__ == "__main__":
    run()
