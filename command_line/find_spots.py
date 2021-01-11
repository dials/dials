# DIALS_ENABLE_COMMAND_LINE_COMPLETION

from __future__ import absolute_import, division, print_function

import logging

from libtbx.phil import parse

from dials.algorithms.shoebox import MaskCode
from dials.algorithms.spot_finding import per_image_analysis
from dials.array_family import flex
from dials.util import log, show_mail_handle_errors
from dials.util.ascii_art import spot_counts_per_image_plot
from dials.util.multi_dataset_handling import generate_experiment_identifiers
from dials.util.options import OptionParser, flatten_experiments
from dials.util.version import dials_version

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
  }

  maximum_trusted_value = None
    .type = float
    .help = "Override maximum trusted value for spot finding only"
    .expert_level = 2

  per_image_statistics = False
    .type = bool
    .help = "Whether or not to print a table of per-image statistics."

  include scope dials.algorithms.spot_finding.factory.phil_scope
""",
    process_includes=True,
)

# Local overrides for dials.find_spots
phil_overrides = parse(
    """
spotfinder {
  mp {
    nproc = Auto
  }
}
"""
)
working_phil = phil_scope.fetch(sources=[phil_overrides])


class Script(object):
    """A class for running the script."""

    def __init__(self, phil=working_phil):
        """Initialise the script."""
        # The script usage
        usage = (
            "usage: dials.find_spots [options] [param.phil] "
            "{models.expt | image1.file [image2.file ...]}"
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

        # Parse the command line
        params, options = self.parser.parse_args(args=args, show_diff_phil=False)

        if __name__ == "__main__":
            # Configure the logging
            log.config(verbosity=options.verbose, logfile=params.output.log)
        logger.info(dials_version())

        # Log the diff phil
        diff_phil = self.parser.diff_phil.as_str()
        if diff_phil != "":
            logger.info("The following parameters have been modified:\n")
            logger.info(diff_phil)

        # Ensure we have a data block
        experiments = flatten_experiments(params.input.experiments)

        # did input have identifier?
        had_identifiers = False
        if all(i != "" for i in experiments.identifiers()):
            had_identifiers = True
        else:
            generate_experiment_identifiers(
                experiments
            )  # add identifier e.g. if coming straight from images

        if len(experiments) == 0:
            self.parser.print_help()
            return

        # If maximum_trusted_value assigned, use this temporarily for the
        # spot finding
        if params.maximum_trusted_value is not None:
            logger.info(
                "Overriding maximum trusted value to %.1f", params.maximum_trusted_value
            )
            input_trusted_ranges = {}
            for _d, detector in enumerate(experiments.detectors()):
                for _p, panel in enumerate(detector):
                    trusted = panel.get_trusted_range()
                    input_trusted_ranges[(_d, _p)] = trusted
                    panel.set_trusted_range((trusted[0], params.maximum_trusted_value))

        # Loop through all the imagesets and find the strong spots
        reflections = flex.reflection_table.from_observations(experiments, params)

        # Add n_signal column - before deleting shoeboxes
        good = MaskCode.Foreground | MaskCode.Valid
        reflections["n_signal"] = reflections["shoebox"].count_mask_values(good)

        # Delete the shoeboxes
        if not params.output.shoeboxes:
            del reflections["shoebox"]

        # ascii spot count per image plot - per imageset

        imagesets = []
        for i, experiment in enumerate(experiments):
            if experiment.imageset not in imagesets:
                imagesets.append(experiment.imageset)

        for imageset in imagesets:
            selected = flex.bool(reflections.nrows(), False)
            for i, experiment in enumerate(experiments):
                if experiment.imageset is not imageset:
                    continue
                selected.set_selected(reflections["id"] == i, True)
            ascii_plot = spot_counts_per_image_plot(reflections.select(selected))
            if len(ascii_plot):
                logger.info("\nHistogram of per-image spot count for imageset %i:" % i)
                logger.info(ascii_plot)

        # Save the reflections to file
        logger.info("\n" + "-" * 80)
        # If started with images and not saving experiments, then remove id mapping
        # as the experiment linked to will no longer exists after exit.
        if not had_identifiers:
            if not params.output.experiments:
                for k in reflections.experiment_identifiers().keys():
                    del reflections.experiment_identifiers()[k]

        reflections.as_file(params.output.reflections)
        logger.info(
            "Saved {} reflections to {}".format(
                len(reflections), params.output.reflections
            )
        )

        # Reset the trusted ranges
        if params.maximum_trusted_value is not None:
            for _d, detector in enumerate(experiments.detectors()):
                for _p, panel in enumerate(detector):
                    trusted = input_trusted_ranges[(_d, _p)]
                    panel.set_trusted_range(trusted)

        # Save the experiments
        if params.output.experiments:

            logger.info("Saving experiments to {}".format(params.output.experiments))
            experiments.as_file(params.output.experiments)

        # Print some per image statistics
        if params.per_image_statistics:
            for i, experiment in enumerate(experiments):
                logger.info("Number of centroids per image for imageset %i:", i)
                refl = reflections.select(reflections["id"] == i)
                refl.centroid_px_to_mm([experiment])
                refl.map_centroids_to_reciprocal_space([experiment])
                stats = per_image_analysis.stats_per_image(
                    experiment, refl, resolution_analysis=False
                )
                logger.info(str(stats))

        if params.output.experiments:
            return experiments, reflections
        else:
            return reflections


@show_mail_handle_errors()
def run(args=None):
    with show_mail_handle_errors():
        script = Script()
        script.run(args)


if __name__ == "__main__":
    run()
