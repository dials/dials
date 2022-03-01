from __future__ import annotations

import logging
import sys

from libtbx.phil import parse

from dials.util import Sorry, log, show_mail_handle_errors

logger = logging.getLogger("dials.command_line.export_best")

help_message = """
This program is used to export the results of dials processing for the strategy
program BEST.
"""

phil_scope = parse(
    """

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

    prefix = best
      .type = str
      .help = "The prefix for the output file names for best"
              "(.hkl, .dat and .par files)"

  }
"""
)


class BestExporter:
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
        prefix = self.params.output.prefix

        best.write_background_file(f"{prefix}.dat", imageset, n_bins=self.params.n_bins)
        best.write_integrated_hkl(prefix, reflections)
        best.write_par_file(f"{prefix}.par", experiment)


@show_mail_handle_errors()
def run(args=None):
    from dials.util.options import (
        ArgumentParser,
        reflections_and_experiments_from_files,
    )
    from dials.util.version import dials_version

    usage = "dials.export models.expt reflections.refl [options]"

    parser = ArgumentParser(
        usage=usage,
        read_experiments=True,
        read_reflections=True,
        phil=phil_scope,
        epilog=help_message,
    )

    # Get the parameters
    params, options = parser.parse_args(args, show_diff_phil=False)

    # Configure the logging
    log.config(logfile=params.output.log)

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
    reflections, experiments = reflections_and_experiments_from_files(
        params.input.reflections, params.input.experiments
    )

    exporter = BestExporter(params, experiments, reflections)
    exporter.export()


if __name__ == "__main__":
    run()
