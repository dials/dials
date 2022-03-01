"""
Estimate a resolution limit based on merging statistics calculated in resolution bins.

A number of metrics are supported for estimating a resolution limit,
including:

- `cc_half` (this is the default)
- `isigma` (unmerged <I/sigI>)
- `misigma` (merged <I/sigI>)
- `i_mean_over_sigma_mean` (unmerged <I>/<sigI>)
- `cc_ref` (CC vs provided reference data set)
- `completeness`
- `rmerge`

Resolution estimation is performed by fitting an appropriate curve to the relevant
merging statistics calculated in resolution bins (with a roughly equal number of
reflections per bin). The estimated resolution limit is chosen as the resolution at
which the fitted function equals the specified criteria.

If multiple metrics are requested, the chosen resolution limit will be the lowest
resolution value estimated across the selected metrics.

The fitting functions for the various metrics are defined as follows:

- `cc_half`: fit a tanh function the form (1/2)(1 - tanh(z)) where z = (s - s0)/r, s0 is
  the value of s at the half-falloff value, and r controls the steepness of falloff
- `isigma`, `misigma`, `i_mean_over_sigma_mean`: fit a polynomial to the values
  log(y(x))
- `rmerge`: fit a polynomial to the values log(1/y(x))
- `completeness`: fit a polynomial to the values y(x)

Example use cases
-----------------

Run with defaults on scaled data::

  dials.estimate_resolution scaled.expt scaled.refl

Run with default on scaled unmerged mtz file::

  dials.estimate_resolution scaled_unmerged.mtz

Override the default cc_half cutoff::

  dials.estimate_resolution scaled.expt scaled.refl cc_half=0.1

Use merged <I/sigI> resolution cutoff instead of cc_half::

  dials.estimate_resolution scaled.expt scaled.refl misigma=1.0 cc_half=None

Use unmerged <I/sigI> resolution cutoff in addition to default cc_half::

  dials.estimate_resolution scaled.expt scaled.refl isigma=0.25

Use cc_ref resolution cutoff::

  dials.estimate_resolution cc_ref=0.3 cc_half=None reference=reference.mtz
"""

from __future__ import annotations

import json
import logging
import sys

from jinja2 import ChoiceLoader, Environment, PackageLoader

import libtbx.phil

from dials.util import log, resolution_analysis, show_mail_handle_errors
from dials.util.multi_dataset_handling import parse_multiple_datasets
from dials.util.options import ArgumentParser, reflections_and_experiments_from_files
from dials.util.version import dials_version

logger = logging.getLogger("dials.estimate_resolution")


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


@show_mail_handle_errors()
def run(args=None):
    usage = "dials.estimate_resolution [options] (scaled.expt scaled.refl | scaled_unmerged.mtz)"

    parser = ArgumentParser(
        usage=usage,
        phil=phil_scope,
        read_reflections=True,
        read_experiments=True,
        check_format=False,
        epilog=__doc__,
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
        if len(experiments) != len(reflections):
            sys.exit(
                f"Mismatched number of experiments and reflection tables found: {len(experiments)} & {len(reflections)}."
            )
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
    run()
