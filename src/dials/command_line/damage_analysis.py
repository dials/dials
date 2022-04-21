"""
The program dials.damage_analysis calculates dose dependent data quality statistics.

The algorithms and statistics are described in "How best to use photons",
G. Winter et al. 2019 https://doi.org/10.1107/S2059798319003528

Input data can be provided in mtz format or scaled .expt and .refl files.

For dials datafiles, by default, the image number is taken as a proxy for dose.
For multi-sweep datasets, each experiment is assumed to be a measurement of a
new crystal, and is assigned a starting dose of zero. If the multi-sweep
dataset is repeated measurements of one crystal, then the option
shared_crystal=True can be given, which will accumulate dose across the sweeps.
Alternatively, a list of starting_doses can be given. The default assumption is
that each image corresponds to an equivalent dose. This can be overwritten with
the dose_per_image option, where a list of values can be provided to set a
different dose per image for each sweep.

For mtz file input, if there is no 'dose' column, then by default the batch
number is used as a proxy for dose (i.e. dose is accumulated across sweeps),
which may not be suitable for multi-sweep mtz files, unless all sweeps are
measured on the same crystal.

Example usage:

dials.damage_analysis scaled.expt scaled.refl

dials.damage_analysis scaled.mtz

dials.damage_analysis scaled.expt scaled.refl shared_crystal=True

"""

from __future__ import annotations

import json
import logging
import os
import sys

from jinja2 import ChoiceLoader, Environment, PackageLoader

from cctbx import crystal, miller
from iotbx import mtz
from libtbx import phil
from scitbx.array_family import flex

from dials.command_line.symmetry import median_unit_cell
from dials.pychef import Statistics, batches_to_dose, interpret_images_to_doses_options
from dials.util import log, resolution_analysis, show_mail_handle_errors
from dials.util.filter_reflections import filter_reflection_table
from dials.util.options import ArgumentParser, reflections_and_experiments_from_files
from dials.util.version import dials_version

try:
    from typing import List
except ImportError:
    pass

logger = logging.getLogger("dials.command_line.damage_analysis")

phil_scope = phil.parse(
    """\
dose {
    experiments
    {
        dose_per_image = 1
        .type = ints(value_min=0)
        .help = "The 'dose' accumulated per image. If more than one value is"
                "given, this indicates the dose per image for each experiment,"
                "and hence must match the number of experiments."
        starting_doses = None
        .type = ints(value_min=0)
        .help = "The dose values at the start of each sweep. Must match the"
                "number of experiments. If none given, it is assumed that each"
                "sweep starts with zero accumulated dose."
        shared_crystal = False
        .type = bool
        .help = "Option to indicate that all sweeps correspond to measurements"
                "on the same crystal. Therefore the starting doses are"
                "automatically adjusted to account for previously accumulated dose."
    }
}
output {
    log = "dials.damage_analysis.log"
        .type = str
        .help = "The log name"
    html = "dials.damage_analysis.html"
        .type = str
        .help = "Filename for the html report."
    json = None
        .type = str
        .help = "Filename for the html report data in json format."
}
include scope dials.pychef.phil_scope
""",
    process_includes=True,
)


class PychefRunner:

    """Class to prepare input data and run the pychef algorithm."""

    def __init__(self, intensities, dose, params):
        self.params = params
        self.stats = None
        sel = dose >= 0
        self.intensities = intensities.select(sel)
        self.dose = dose.select(sel)
        self.resolution_filter()

    def resolution_filter(self):
        """Filter arrays on resolution."""
        if not self.params.d_min and self.params.min_completeness:
            # Update self.params.d_min using resolution estimate
            params = resolution_analysis.phil_defaults.extract().resolution
            params.nbins = self.params.resolution_bins
            r = resolution_analysis.Resolutionizer(self.intensities, params)
            self.params.d_min = r.resolution(
                resolution_analysis.metrics.COMPLETENESS,
                limit=self.params.min_completeness,
            ).d_min
            logger.info("Estimated d_min: %.2f", self.params.d_min)

        if self.params.d_min or self.params.d_max:
            sel = flex.bool(self.intensities.size(), True)
            d_spacings = self.intensities.d_spacings().data()
            if self.params.d_min:
                sel &= d_spacings >= self.params.d_min
            if self.params.d_max:
                sel &= d_spacings <= self.params.d_max
            self.intensities = self.intensities.select(sel)
            self.dose = self.dose.select(sel)

    @classmethod
    def from_mtz(cls, params, mtz_object):
        """Initialise the class from mtz input.

        Args:
            params: A damage-analysis phil params object
            mtz_object: An iotbx.mtz.object.
        """
        miller_arrays = mtz_object.as_miller_arrays(
            merge_equivalents=False, anomalous=params.anomalous
        )

        intensities = None
        batches = None
        dose = None

        for ma in miller_arrays:
            if ma.info().labels == ["BATCH"]:
                batches = ma
            elif ma.info().labels == ["DOSE"]:
                dose = ma
            elif ma.info().labels == ["I", "SIGI"]:
                intensities = ma
            elif ma.info().labels == ["I(+)", "SIGI(+)", "I(-)", "SIGI(-)"]:
                intensities = ma
        if not intensities:
            raise KeyError("Intensity array not found in mtz file")
        if not batches:
            raise KeyError("Batch array not found in mtz file")

        indices = mtz_object.extract_original_index_miller_indices()
        intensities = intensities.customized_copy(indices=indices)
        batches = batches.customized_copy(indices=indices)

        if params.anomalous:
            intensities = intensities.as_anomalous_array()
            batches = batches.as_anomalous_array()

        if dose is None:
            dose = batches_to_dose(batches.data(), params.dose)
        else:
            dose = dose.data()

        return cls(intensities, dose, params)

    @classmethod
    def from_dials_data_files(cls, params, experiments, reflection_table):
        """Initialise the class from an experiment list and reflection table.

        Args:
            params: A damage-analysis phil params object
            experiments: An ExperimentList
            reflection_table: A reflection table.
        """
        reflection_table = filter_reflection_table(
            reflection_table, intensity_choice=["scale"], partiality_threshold=0.4
        )

        # get scaled intensities
        intensities = miller.array(
            miller.set(
                crystal.symmetry(
                    unit_cell=median_unit_cell(experiments),
                    space_group=experiments[0].crystal.get_space_group(),
                ),
                indices=reflection_table["miller_index"],
                anomalous_flag=params.anomalous,
            ),
            data=reflection_table["intensity.scale.value"],
            sigmas=flex.sqrt(reflection_table["intensity.scale.variance"]),
        )
        intensities.set_observation_type_xray_intensity()

        doses = flex.double()
        start_doses, doses_per_image = interpret_images_to_doses_options(
            experiments,
            params.dose.experiments.dose_per_image,
            params.dose.experiments.starting_doses,
            params.dose.experiments.shared_crystal,
        )
        logger.info(
            "Interpreting data using:\n  starting_doses=%s\n  dose_per_image=%s",
            ", ".join("%s" % i for i in start_doses)
            if len(set(start_doses)) > 1
            else f" all {start_doses[0]}",
            ", ".join("%s" % i for i in doses_per_image)
            if len(set(doses_per_image)) > 1
            else f" all {doses_per_image[0]}",
        )

        for expt, starting_dose, dose_per_img in zip(
            experiments, start_doses, doses_per_image
        ):
            refls = reflection_table.select(expt)
            imgno = flex.ceil(refls["xyzobs.px.value"].parts()[2])
            dose = (imgno * dose_per_img) + starting_dose
            doses.extend(dose)

        doses = doses.iround()

        return cls(intensities, doses, params)

    def run(self):
        """Run the pychef analysis."""
        self.stats = Statistics(
            self.intensities,
            self.dose,
            n_bins=self.params.resolution_bins,
            range_min=self.params.range.min,
            range_max=self.params.range.max,
            range_width=self.params.range.width,
        )

        logger.debug(self.stats.completeness_vs_dose_str())
        logger.debug(self.stats.rcp_vs_dose_str())
        logger.debug(self.stats.scp_vs_dose_str())
        logger.debug(self.stats.rd_vs_dose_str())

    def make_html_report(self, html_filename=None, json_filename=None):
        """Generate html report from pychef stats."""
        data = {"dose_plots": self.stats.to_dict()}
        if html_filename:
            logger.info("Writing html report to: %s", html_filename)
            loader = ChoiceLoader(
                [
                    PackageLoader("dials", "templates"),
                    PackageLoader("dials", "static", encoding="utf-8"),
                ]
            )
            env = Environment(loader=loader)
            template = env.get_template("simple_report.html")
            html = template.render(
                page_title="Damage analysis report",
                panel_title="Damage analysis plots",
                panel_id="dose_plots",
                graphs=data["dose_plots"],
            )
            with open(html_filename, "wb") as f:
                f.write(html.encode("utf-8", "xmlcharrefreplace"))
        if json_filename:
            logger.info("Writing html report data to: %s", json_filename)
            with open(json_filename, "w") as outfile:
                json.dump(data, outfile)


@show_mail_handle_errors()
def run(args: List[str] = None, phil: phil.scope = phil_scope) -> None:
    """Run the command-line script."""

    usage = "dials.damage_analysis [options] scaled.expt scaled.refl | scaled.mtz"

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

    reflections, experiments = reflections_and_experiments_from_files(
        params.input.reflections, params.input.experiments
    )
    try:
        if experiments and reflections:
            if len(reflections) != 1:
                raise ValueError("A single input reflections datafile is required")
            if "inverse_scale_factor" not in reflections[0]:
                raise KeyError("Input data must be scaled.")
            script = PychefRunner.from_dials_data_files(
                params,
                experiments,
                reflections[0],
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
                script = PychefRunner.from_mtz(params, mtz_object)
        else:
            parser.print_help()
            raise ValueError("Suitable input datafiles not provided")
    except (ValueError, KeyError) as e:
        sys.exit(f"Error: {e}")
    else:
        script.run()
        script.make_html_report(params.output.html, params.output.json)


if __name__ == "__main__":
    run()
