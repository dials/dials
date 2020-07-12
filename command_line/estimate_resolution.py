# LIBTBX_SET_DISPATCHER_NAME dials.estimate_resolution
# LIBTBX_SET_DISPATCHER_NAME dials.resolutionizer
from __future__ import absolute_import, division, print_function

import json
import logging
import sys
import warnings
from jinja2 import Environment, ChoiceLoader, PackageLoader

import libtbx.phil

from dials.util import resolution_analysis
from dials.util import log
from dials.util.options import OptionParser, reflections_and_experiments_from_files
from dials.util.version import dials_version
from dials.util.multi_dataset_handling import parse_multiple_datasets

logger = logging.getLogger("dials.estimate_resolution")

help_message = """
"""


phil_scope = libtbx.phil.parse(
    """
include scope dials.util.resolution_analysis.phil_defaults

output {
  log = dials.estimate_resolution.log
    .type = path
  html = dials.estimate_resolution.html
    .type = path
  json = None
    .type = path
}
""",
    process_includes=True,
)


def run(args):
    usage = "dials.estimate_resolution [options] (scaled.expt scaled.refl | scaled_unmerged.mtz)"

    import libtbx.load_env

    if libtbx.env.dispatcher_name == "dials.resolutionizer":
        warnings.warn(
            "dials.resolutionizer is now deprecated, please use dials.estimate_resolution instead"
        )

    parser = OptionParser(
        usage=usage,
        phil=phil_scope,
        read_reflections=True,
        read_experiments=True,
        check_format=False,
        epilog=help_message,
    )

    params, options, unhandled = parser.parse_args(
        args=args, return_unhandled=True, show_diff_phil=True
    )

    reflections, experiments = reflections_and_experiments_from_files(
        params.input.reflections, params.input.experiments
    )
    if (not reflections or not experiments) and not unhandled:
        parser.print_help()
        return

    if reflections and experiments and unhandled:
        sys.exit(
            "Must provide either scaled unmerged mtz OR dials-format scaled reflections and experiments files"
        )

    # Configure the logging
    log.config(logfile=params.output.log, verbosity=options.verbose)
    logger.info(dials_version())

    if len(unhandled) == 1:
        scaled_unmerged = unhandled[0]
        m = resolution_analysis.Resolutionizer.from_unmerged_mtz(
            scaled_unmerged, params.resolution
        )
    else:
        reflections = parse_multiple_datasets(reflections)
        m = resolution_analysis.Resolutionizer.from_reflections_and_experiments(
            reflections, experiments, params.resolution
        )

    plots = m.resolution_auto()

    if params.output.html:
        output_html_report(plots, params.output.html)

    if params.output.json:
        with open(params.output.json, "w") as fh:
            json.dump(plots, fh)

    return plots


def output_html_report(plots, filename):
    loader = ChoiceLoader(
        [
            PackageLoader("dials", "templates"),
            PackageLoader("dials", "static", encoding="utf-8"),
        ]
    )
    env = Environment(loader=loader)
    template = env.get_template("simple_report.html")
    html = template.render(
        page_title="dials.estimate_resolution report",
        panel_title="Analysis by resolution",
        panel_id="dials_estimate_resolution",
        graphs=plots,
    )
    with open(filename, "wb") as f:
        f.write(html.encode("utf-8", "xmlcharrefreplace"))


if __name__ == "__main__":
    run(sys.argv[1:])
