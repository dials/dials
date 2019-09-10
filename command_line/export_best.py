from __future__ import absolute_import, division, print_function

import logging
import sys

from libtbx.phil import parse
from libtbx import Auto

logger = logging.getLogger("dials.command_line.export_best")

help_message = """
This program is used to export the results of dials processing for the strategy
program BEST.

"""

phil_scope = parse(
    """

  intensity = *auto profile sum scale
    .type = choice(multi=True)
    .help = "Choice of which intensities to export. Allowed combinations:
            scale, profile, sum, profile+sum, sum+profile+scale. Auto will
            default to scale or profile+sum depending on if the data are scaled."

  debug = False
    .type = bool
    .help = "Output additional debugging information"

  prefix = best
    .type = str
    .help = "The prefix for the output file names for best"
            "(.hkl, .dat and .par files)"

  n_bins = 100
    .type = int(value_min=1)
    .help = "Number of resolution bins for background estimation"

  min_partiality = 0.1
    .type = float(value_min=0, value_max=1)
    .help = "Minimum partiality of reflections to export"

  output {

    log = dials.export_best.log
      .type = path
      .help = "The log filename"

  }
"""
)


class BestExporter(object):
    def __init__(self, params, experiments, reflections):
        """
        Initialise the exporter

        :param params: The phil parameters
        :param experiments: The experiment list
        :param reflections: The reflection tables
        """

        # Check the input
        if not experiments:
            raise Sorry("BEST exporter requires an experiment list")
        if not reflections:
            raise Sorry("BEST exporter require a reflection table")

        # Save the stuff
        self.params = params
        self.experiments = experiments
        self.reflections = reflections

    def export(self):
        """
        Export the files
        """
        from dials.util import best

        experiment = self.experiments[0]
        reflections = self.reflections[0]
        partiality = reflections["partiality"]
        sel = partiality >= self.params.min_partiality
        logger.info(
            "Selecting %s/%s reflections with partiality >= %s",
            sel.count(True),
            sel.size(),
            self.params.min_partiality,
        )
        if sel.count(True) == 0:
            raise Sorry(
                "No reflections remaining after filtering for minimum partiality (min_partiality=%f)"
                % (self.params.min_partiality)
            )
        reflections = reflections.select(sel)

        imageset = experiment.imageset
        prefix = self.params.prefix

        best.write_background_file(
            "%s.dat" % prefix, imageset, n_bins=self.params.n_bins
        )
        best.write_integrated_hkl(prefix, reflections)
        best.write_par_file("%s.par" % prefix, experiment)


if __name__ == "__main__":
    from dials.util.options import (
        OptionParser,
        flatten_experiments,
        flatten_reflections,
    )
    from dials.util.version import dials_version
    from dials.util import log
    from dials.util import Sorry

    usage = "dials.export models.expt reflections.pickle [options]"

    parser = OptionParser(
        usage=usage,
        read_experiments=True,
        read_reflections=True,
        phil=phil_scope,
        epilog=help_message,
    )

    # Get the parameters
    params, options = parser.parse_args(show_diff_phil=False)

    # Configure the logging
    log.config(info=params.output.log)

    # Print the version number
    logger.info(dials_version())

    # Log the diff phil
    diff_phil = parser.diff_phil.as_str()
    if diff_phil != "":
        logger.info("The following parameters have been modified:\n")
        logger.info(diff_phil)

    if not params.input.experiments and not params.input.reflections:
        parser.print_help()
        sys.exit()

    # Get the experiments and reflections
    experiments = flatten_experiments(params.input.experiments)
    reflections = flatten_reflections(params.input.reflections)

    # do auto intepreting of intensity choice:
    # note that this may still fail certain checks further down the processing,
    # but these are the defaults to try
    if params.intensity in ([None], [Auto], ["auto"]) and reflections:
        if ("intensity.scale.value" in reflections[0]) and (
            "intensity.scale.variance" in reflections[0]
        ):
            params.intensity = ["scale"]
            logger.info("Data appears to be scaled, setting intensity = scale")
        else:
            params.intensity = ["profile", "sum"]
            logger.info("Data appears to be unscaled, setting intensity = profile+sum")

    exporter = BestExporter(params, experiments, reflections)
    exporter.export()
