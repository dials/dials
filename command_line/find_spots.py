#!/usr/bin/env python
#
# dials.find_spots.py
#
#  Copyright (C) 2013 Diamond Light Source
#
#  Author: James Parkhurst
#
#  This code is distributed under the BSD license, a copy of which is
#  included in the root directory of this package.
#
# DIALS_ENABLE_COMMAND_LINE_COMPLETION

from __future__ import absolute_import, division, print_function

import logging

logger = logging.getLogger("dials.command_line.find_spots")

help_message = """

This program tries to find strong spots on a sequence of images. The program can
be called with either a "models.expt" file or a sequence of image files (see
help for dials.import for more information about how images are imported). Spot
finding will be done against each logically grouped set of images given. Strong
pixels will be found on each image and spots will be formed from connected
components. In the case of rotation images, connected component labelling will
be done in 3D.

Once a set of spots have been found, their centroids and intensities will be
calculated. They will then be filtered according to the particular preferences
of the user. The output will be a file (strong.refl) containing a list of spot
centroids and intensities which can be used in the dials.index program. To view
a list of parameters for spot finding use the --show-config option.

Examples::

  dials.find_spots image1.cbf

  dials.find_spots imager_00*.cbf

  dials.find_spots models.expt

  dials.find_spots models.expt output.reflections=strong.refl

"""

# Set the phil scope
from libtbx.phil import parse

phil_scope = parse(
    """

  output {
    reflections = 'strong.refl'
      .type = str
      .help = "The output filename"

    shoeboxes = True
      .type = bool
      .help = "Save the raw pixel values inside the reflection shoeboxes."

    experiments = None
      .type = str
      .help = "Save the modified experiments."
              "(usually only modified with hot pixel mask)"

    log = 'dials.find_spots.log'
      .type = str
      .help = "The log filename"

    debug_log = 'dials.find_spots.debug.log'
      .type = str
      .help = "The debug log filename"
  }

  per_image_statistics = False
    .type = bool
    .help = "Whether or not to print a table of per-image statistics."

  verbosity = 0
    .type = int(value_min=0)
    .help = "The verbosity level"

  include scope dials.algorithms.spot_finding.factory.phil_scope

""",
    process_includes=True,
)


class Script(object):
    """A class for running the script."""

    def __init__(self, phil=phil_scope):
        """Initialise the script."""
        from dials.util.options import OptionParser
        import libtbx.load_env

        # The script usage
        usage = (
            "usage: %s [options] [param.phil] "
            "{models.expt | image1.file [image2.file ...]}" % libtbx.env.dispatcher_name
        )

        # Initialise the base class
        self.parser = OptionParser(
            usage=usage,
            phil=phil,
            epilog=help_message,
            read_experiments_from_images=True,
            read_experiments=True,
        )

    def run(self, args=None):
        """Execute the script."""
        from dxtbx.model.experiment_list import ExperimentListDumper
        from dials.array_family import flex
        from dials.util.options import flatten_experiments
        from time import time
        from dials.util import log

        start_time = time()

        # Parse the command line
        params, options = self.parser.parse_args(args=args, show_diff_phil=False)

        if __name__ == "__main__":
            # Configure the logging
            log.config(
                params.verbosity, info=params.output.log, debug=params.output.debug_log
            )

        from dials.util.version import dials_version

        logger.info(dials_version())

        # Log the diff phil
        diff_phil = self.parser.diff_phil.as_str()
        if diff_phil != "":
            logger.info("The following parameters have been modified:\n")
            logger.info(diff_phil)

        # Ensure we have a data block
        experiments = flatten_experiments(params.input.experiments)
        if len(experiments) == 0:
            self.parser.print_help()
            return

        # Loop through all the imagesets and find the strong spots
        reflections = flex.reflection_table.from_observations(experiments, params)

        # Add n_signal column - before deleting shoeboxes
        from dials.algorithms.shoebox import MaskCode

        good = MaskCode.Foreground | MaskCode.Valid
        reflections["n_signal"] = reflections["shoebox"].count_mask_values(good)

        # Delete the shoeboxes
        if not params.output.shoeboxes:
            del reflections["shoebox"]

        # ascii spot count per image plot
        from dials.util.ascii_art import spot_counts_per_image_plot

        for i, experiment in enumerate(experiments):
            ascii_plot = spot_counts_per_image_plot(
                reflections.select(reflections["id"] == i)
            )
            if len(ascii_plot):
                logger.info("\nHistogram of per-image spot count for imageset %i:" % i)
                logger.info(ascii_plot)

        # Save the reflections to file
        logger.info("\n" + "-" * 80)
        reflections.as_pickle(params.output.reflections)
        logger.info(
            "Saved {} reflections to {}".format(
                len(reflections), params.output.reflections
            )
        )

        # Save the experiments
        if params.output.experiments:
            logger.info("Saving experiments to {}".format(params.output.experiments))
            dump = ExperimentListDumper(experiments)
            dump.as_file(params.output.experiments)

        # Print some per image statistics
        if params.per_image_statistics:
            from dials.algorithms.spot_finding import per_image_analysis
            from six.moves import cStringIO as StringIO

            s = StringIO()
            for i, experiment in enumerate(experiments):
                print("Number of centroids per image for imageset %i:" % i, file=s)
                imageset = experiment.imageset
                stats = per_image_analysis.stats_imageset(
                    imageset,
                    reflections.select(reflections["id"] == i),
                    resolution_analysis=False,
                )
                per_image_analysis.print_table(stats, out=s)
            logger.info(s.getvalue())

        # Print the time
        logger.info("Time Taken: %f" % (time() - start_time))

        if params.output.experiments:
            return experiments, reflections
        else:
            return reflections


if __name__ == "__main__":
    from dials.util import halraiser

    try:
        script = Script()
        script.run()
    except Exception as e:
        halraiser(e)
