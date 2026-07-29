from __future__ import annotations

import json
import logging
import os
import sys

from jinja2 import ChoiceLoader, Environment, PackageLoader

from iotbx import mtz
from libtbx import phil

from dials.algorithms.scaling.scaling_library import (
    merging_stats_from_scaled_array,
    scaled_data_as_miller_array,
)
from dials.util import log, show_mail_handle_errors
from dials.util.options import ArgumentParser, reflections_and_experiments_from_files
from dials.util.version import dials_version

try:
    from typing import List
except ImportError:
    pass
import iotbx.phil

from dials.report.analysis import (
    table_1_summary,
)

phil_scope = iotbx.phil.parse("""\
n_bins = 10
    .type = int
    .help = "The number of resolution bins for merging statistics calculations"
calculate_sigma_tau_stats = False
    .type = bool
    .help = "Additionally calculate CC1/2 statistics with the sigma-tau method."
output.log=dials.merging_statistics.log
    .type = str
    .help = "The log filename"
output.json=dials.merging_statistics.json
    .type = str
    .help = "The json filename for the plot data."
""")


logger = logging.getLogger("dials")


def generate_html_report(json_data, filename):
    multi_data = None
    styles = {}
    if "multi_data" in json_data:
        multi_data = json_data.pop("multi_data")
    for wl, v in json_data.items():
        str_wl = wl.replace(".", "_")
        for plot_cat in [
            "resolution_plots",
            "misc_plots",
            "unit_cell_plots",
            "orientation_graphs",
        ]:
            if plot_cat in v:
                for name in list(v[plot_cat].keys()):
                    v[plot_cat][name + "_" + str_wl] = v[plot_cat].pop(name)
                    if plot_cat == "orientation_graphs":
                        styles[name + "_" + str_wl] = "square-plot"
            else:
                v[plot_cat] = {}

    loader = ChoiceLoader(
        [
            PackageLoader("dials", "templates"),
            PackageLoader("dials", "static", encoding="utf-8"),
        ]
    )
    env = Environment(loader=loader)
    template = env.get_template("merging_statistics_report.html")
    html = template.render(
        page_title="dials merging statistics report",
        individual_reports=json_data,
        multi_data=multi_data,
        styles=styles,
    )
    logger.info("Writing html report to %s", filename)
    with open(filename, "wb") as f:
        f.write(html.encode("utf-8", "xmlcharrefreplace"))


@show_mail_handle_errors()
def run(args: List[str] = None, phil: phil.scope = phil_scope) -> None:
    """Run the command-line script."""

    usage = "dials.merging_statistics [options] scaled.expt scaled.refl | scaled.mtz"

    parser = ArgumentParser(
        usage=usage,
        phil=phil,
        epilog=__doc__,
        read_experiments=True,
        read_reflections=True,
        check_format=False,
    )

    params, _, unhandled = parser.parse_args(
        args=args, show_diff_phil=False, return_unhandled=True
    )

    log.config(logfile=params.output.log)
    logger.info(dials_version())

    diff_phil = parser.diff_phil.as_str()
    if diff_phil:
        logger.info("The following parameters have been modified:\n%s", diff_phil)

    reflections, experiments = reflections_and_experiments_from_files(
        params.input.reflections, params.input.experiments
    )
    try:
        if experiments and reflections:
            if len(reflections) != 1:
                raise ValueError("A single input reflections datafile is required")
            if "inverse_scale_factor" not in reflections[0]:
                raise KeyError("Input data must be scaled.")
            iobs = scaled_data_as_miller_array(
                reflections, experiments, anomalous_flag=True
            )

        elif unhandled and os.path.isfile(unhandled[0]):
            try:
                mtz_object = mtz.object(file_name=unhandled[0])
            except RuntimeError:
                # If an error is encountered trying to read the file as an mtzfile
                raise ValueError(
                    "Input file cannot be read as a valid experiment/reflection file or MTZ file"
                )
            else:
                miller_arrays = mtz_object.as_miller_arrays(
                    merge_equivalents=False, anomalous=True
                )
                intensities = None
                for ma in miller_arrays:
                    if ma.info().labels == ["I", "SIGI"]:
                        intensities = ma
                    elif ma.info().labels == ["I", "SigI"]:
                        intensities = ma
                    elif ma.info().labels == ["IOBS", "SIGIOBS"]:
                        intensities = ma
                    elif ma.info().labels == ["I(+)", "SIGI(+)", "I(-)", "SIGI(-)"]:
                        intensities = ma
                if not intensities:
                    raise KeyError("Intensity array not found in mtz file")
                indices = mtz_object.extract_original_index_miller_indices()
                iobs = intensities.customized_copy(
                    indices=indices, info=intensities.info(), anomalous_flag=True
                )
    except (ValueError, KeyError) as e:
        sys.exit(f"Error: {e}")
    else:
        stats, anom_stats = merging_stats_from_scaled_array(
            iobs,
            additional_stats=True,
            n_bins=params.n_bins,
            sigma_tau_stats=params.calculate_sigma_tau_stats,
        )
        from dials.report.plots import ResolutionPlotsAndStats

        logger.info(table_1_summary(stats, anom_stats))

        is_centric = iobs.space_group().is_centric()
        plotter = ResolutionPlotsAndStats(stats, anom_stats, is_centric)
        d = {
            "scaling_tables": {},
            "resolution_plots": {},
            "batch_plots": {},
            "misc_plots": {},
            "anom_plots": {},
            "image_range_tables": [],
        }
        d["resolution_plots"].update(plotter.make_all_plots())
        json_data = {"main": d}
        generate_html_report(json_data, "dials.merging_statistics.html")
        with open(params.output.json, "w") as f:
            json.dump(json_data, f, indent=2)


if __name__ == "__main__":
    run()
