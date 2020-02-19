# coding: utf-8
"""
The program dials.dose_analysis calculates dose dependent data quality statistics.

The algorithms and statistics are described in "How best to use photons",
G. Winter et al. 2019 https://doi.org/10.1107/S2059798319003528

The program makes a noteable assumption which may be unsuitable for particular
experiments.
For running on dials datafiles (expt, refl), the image number is
taken as a proxy for dose (calculated individually on each sweep for multi-
sweep datasets).
For mtz file input, if there is no 'dose' column, then the batch number is
used as a proxy for dose, which may not be suitable for multi-sweep mtz files,
unless all sweeps are measured on the same crystal.

"""
from __future__ import absolute_import, division, print_function
import logging
import sys
from libtbx import phil
from dials.util import log, show_mail_on_error
from dials.util.options import OptionParser, reflections_and_experiments_from_files
from dials.util.version import dials_version
from dials.util.filter_reflections import filter_reflection_table
from dials.util.resolutionizer import Resolutionizer, phil_defaults
from dials.command_line.symmetry import median_unit_cell
from dials.pychef import batches_to_dose, Statistics
from iotbx import mtz
from scitbx.array_family import flex
from cctbx import miller, crystal
from jinja2 import Environment, ChoiceLoader, PackageLoader

try:
    from typing import List
except ImportError:
    pass

logger = logging.getLogger("dials.command_line.dose_analysis")

phil_scope = phil.parse(
    """\
input {
    mtzfile = None
        .type = str
        .help = "We can also import an MTZ file"
}
output {
    log = "dials.dose_analysis.log"
        .type = str
        .help = "The log name"
    html = "dials.dose_analysis.html"
        .type = str
        .help = "Filename for the html report."
}
include scope dials.pychef.phil_scope
""",
    process_includes=True,
)


class PychefRunner(object):

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
            params = phil_defaults.extract().resolutionizer
            params.nbins = self.params.resolution_bins
            r = Resolutionizer(self.intensities, params)
            self.params.d_min = r.resolution_completeness(
                limit=self.params.min_completeness
            )
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
    def from_mtzfile(cls, params, mtzfile):
        """Initialise the class from mtzfile input.

        Args:
            params: A dose-analysis phil params object
            mtzfile: Filename of the mtz file.
        """
        mtz_object = mtz.object(file_name=mtzfile)
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
    def from_dials_datafiles(cls, params, experiments, reflection_table):
        """Initialise the class from an experiment list and reflection table.

        Args:
            params: A dose-analysis phil params object
            experiments: An ExperimentList
            reflection_table: A reflection table.
        """
        reflection_table = filter_reflection_table(
            reflection_table, intensity_choice=["scale"]
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
            sigmas=reflection_table["intensity.scale.variance"] ** 0.5,
        )
        intensities.set_observation_type_xray_intensity()
        dose = reflection_table["xyzobs.px.value"].parts()[2]
        dose = dose.iround()

        return cls(intensities, dose, params)

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

    def make_html_report(self):
        """Generate html report from pychef stats."""
        data = {"dose_plots": self.stats.to_dict()}

        logger.info("Writing html report to: %s", self.params.output.html)
        loader = ChoiceLoader(
            [
                PackageLoader("dials", "templates"),
                PackageLoader("dials", "static", encoding="utf-8"),
            ]
        )
        env = Environment(loader=loader)
        template = env.get_template("dose_analysis_report.html")
        html = template.render(
            page_title="Dose analysis report", dose_plots=data["dose_plots"]
        )
        with open(self.params.output.html, "wb") as f:
            f.write(html.encode("utf-8", "xmlcharrefreplace"))


def run(args=None, phil=phil_scope):  # type: (List[str], phil.scope) -> None
    """Run the command-line script."""

    usage = "dials.dose_analysis [options] scaled.expt scaled.refl | mtzfile=scaled.mtz"

    parser = OptionParser(
        usage=usage,
        phil=phil,
        epilog=__doc__,
        read_experiments=True,
        read_reflections=True,
        check_format=False,
    )

    params, _ = parser.parse_args(args=args, show_diff_phil=False)

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
            script = PychefRunner.from_dials_datafiles(
                params, experiments, reflections[0],
            )
        elif params.input.mtzfile:
            script = PychefRunner.from_mtzfile(params, params.input.mtzfile)
        else:
            parser.print_help()
            raise ValueError("Suitable input datafiles not provided")
    except (ValueError, KeyError) as e:
        sys.exit("Error: %s" % e.message)
    else:
        script.run()
        script.make_html_report()


if __name__ == "__main__":
    with show_mail_on_error():
        run()
