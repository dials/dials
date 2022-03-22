"""
A standalone program to explore error model minimisation after scaling.

Scale factors and outliers determined from scaling are respected.
Profile intensities will be used by default. Profile/summation intensity
combination can be used, with the options:
    intensity_choice=combine combine.Imid={value}

All scaling corrections are applied to the intensities and variances (except
of course the error model adjustment) before the analysis is done.
"""

from __future__ import annotations

import json
import logging
import sys

from jinja2 import ChoiceLoader, Environment, PackageLoader

import libtbx.phil
from dxtbx import flumpy
from iotbx import phil

from dials.algorithms.scaling.combine_intensities import combine_intensities
from dials.algorithms.scaling.error_model.engine import run_error_model_refinement
from dials.algorithms.scaling.error_model.error_model import (
    BasicErrorModel,
    calc_deltahl,
    calc_sigmaprime,
)
from dials.algorithms.scaling.error_model.error_model_target import (
    calculate_regression_x_y,
)
from dials.algorithms.scaling.Ih_table import IhTable
from dials.algorithms.scaling.plots import (
    error_model_variance_plot,
    error_regression_plot,
    normal_probability_plot,
)
from dials.algorithms.scaling.scaling_library import choose_initial_scaling_intensities
from dials.algorithms.scaling.scaling_utilities import calculate_prescaling_correction
from dials.report.plots import i_over_sig_i_vs_i_plot
from dials.util import log, show_mail_handle_errors
from dials.util.options import ArgumentParser, reflections_and_experiments_from_files
from dials.util.version import dials_version

try:
    from typing import List
except ImportError:
    pass

logger = logging.getLogger("dials")

phil_scope = phil.parse(
    """
    include scope dials.algorithms.scaling.error_model.error_model.phil_scope
    intensity_choice = *profile sum combine
        .type = choice
        .help = "Use profile or summation intensities"
    combine.Imid = None
        .type = float
        .help = "Midpoint value to use when combining profile/summation intensities"
    output {
        log = dials.refine_error_model.log
            .type = str
            .help = "The log filename"
        html = error_model.html
            .type = str
            .help = "Filename for html report of error model"
        json = None
            .type = str
            .help = "Filename for json export of html report data"
    }
    """,
    process_includes=True,
)


def refine_error_model(params, experiments, reflection_tables):
    """Do error model refinement."""

    # prepare relevant data for datastructures
    for i, table in enumerate(reflection_tables):
        # First get the good data
        table = table.select(~table.get_flags(table.flags.bad_for_scaling, all=False))

        # Now chose intensities, ideally these two options could be combined
        # with a smart refactor
        if params.intensity_choice == "combine":
            if not params.combine.Imid:
                sys.exit("Imid value must be provided if intensity_choice=combine")
            table = calculate_prescaling_correction(table)  # needed for below.
            I, V = combine_intensities(table, params.combine.Imid)
            table["intensity"] = I
            table["variance"] = V
        else:
            table = choose_initial_scaling_intensities(
                table, intensity_choice=params.intensity_choice
            )
        reflection_tables[i] = table
    space_group = experiments[0].crystal.get_space_group()
    Ih_table = IhTable(
        reflection_tables, space_group, additional_cols=["partiality"], anomalous=True
    )

    # now do the error model refinement
    model = BasicErrorModel(basic_params=params.basic)
    try:
        model = run_error_model_refinement(model, Ih_table)
    except (ValueError, RuntimeError) as e:
        logger.info(e)
    else:
        return model


def make_output(model, params):
    """Get relevant data from the model and make html report."""

    if not (params.output.html or params.output.json):
        return

    data = {}
    table = model.filtered_Ih_table
    data["intensity"] = flumpy.from_numpy(table.intensities)
    sigmaprime = calc_sigmaprime(model.parameters, table)
    data["delta_hl"] = calc_deltahl(table, table.calc_nh(), sigmaprime)
    data["inv_scale"] = table.inverse_scale_factors
    data["sigma"] = flumpy.from_numpy(sigmaprime * table.inverse_scale_factors)
    data["binning_info"] = model.binner.binning_info
    d = {"error_model_plots": {}}
    d["error_model_plots"].update(normal_probability_plot(data))
    d["error_model_plots"].update(
        i_over_sig_i_vs_i_plot(data["intensity"], data["sigma"])
    )
    d["error_model_plots"].update(error_model_variance_plot(data))

    if params.basic.minimisation == "regression":
        x, y = calculate_regression_x_y(model.filtered_Ih_table)
        data["regression_x"] = x
        data["regression_y"] = y
        data["model_a"] = model.parameters[0]
        data["model_b"] = model.parameters[1]
        d["error_model_plots"].update(error_regression_plot(data))

    if params.output.html:
        logger.info("Writing html report to: %s", params.output.html)
        loader = ChoiceLoader(
            [
                PackageLoader("dials", "templates"),
                PackageLoader("dials", "static", encoding="utf-8"),
            ]
        )
        env = Environment(loader=loader)
        template = env.get_template("simple_report.html")
        html = template.render(
            page_title="DIALS error model refinement report",
            panel_title="Error distribution plots",
            panel_id="error_model_plots",
            graphs=d["error_model_plots"],
        )
        with open(params.output.html, "wb") as f:
            f.write(html.encode("utf-8", "xmlcharrefreplace"))

    if params.output.json:
        logger.info("Writing html report data to: %s", params.output.json)
        with open(params.output.json, "w") as outfile:
            json.dump(d, outfile)


@show_mail_handle_errors()
def run(args: List[str] = None, phil: libtbx.phil.scope = phil_scope) -> None:
    """Run the scaling from the command-line."""
    usage = """Usage: dials.refine_error_model scaled.refl scaled.expt [options]"""

    parser = ArgumentParser(
        usage=usage,
        read_experiments=True,
        read_reflections=True,
        phil=phil,
        check_format=False,
        epilog=__doc__,
    )
    params, options = parser.parse_args(args=args, show_diff_phil=False)

    if not params.input.experiments or not params.input.reflections:
        parser.print_help()
        sys.exit()

    reflections, experiments = reflections_and_experiments_from_files(
        params.input.reflections, params.input.experiments
    )

    log.config(verbosity=options.verbose, logfile=params.output.log)
    logger.info(dials_version())

    diff_phil = parser.diff_phil.as_str()
    if diff_phil:
        logger.info("The following parameters have been modified:\n%s", diff_phil)

    model = refine_error_model(params, experiments, reflections)

    if model:
        make_output(model, params)
    logger.info("Finished")


if __name__ == "__main__":
    run()
