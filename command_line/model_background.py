#!/usr/bin/env python
#
# dials.model_background.py
#
#  Copyright (C) 2013 Diamond Light Source
#
#  Author: James Parkhurst
#
#  This code is distributed under the BSD license, a copy of which is
#  included in the root directory of this package.
#

from __future__ import absolute_import, division, print_function

import logging

logger = logging.getLogger("dials.command_line.model_background")

help_message = """


"""

# Set the phil scope
from libtbx.phil import parse

phil_scope = parse(
    """

  output {
    model = 'background.pickle'
      .type = str
      .help = "The output filename"

    log = 'dials.model_background.log'
      .type = str
      .help = "The log filename"

    debug_log = 'dials.model_background.debug.log'
      .type = str
      .help = "The debug log filename"

    mean_image_prefix = 'mean'
      .type = str
      .help = "The mean background image"

    variance_image_prefix = 'variance'
      .type = str
      .help = "The variance background image"

    dispersion_image_prefix = 'dispersion'
      .type = str
      .help = "The dispersion background image"

    mask_image_prefix = 'mask'
      .type = str
      .help = "The mask background image"

    min_image_prefix = 'min'
      .type = str
      .help = "The min background image"

    max_image_prefix = 'max'
      .type = str
      .help = "The max background image"

    model_image_prefix = 'model'
      .type = str
      .help = "The model background image"

    polar_model_image_prefix = 'polar'
      .type = str
      .help = "The polar model background image"
  }

  verbosity = 0
    .type = int(value_min=0)
    .help = "The verbosity level"

  modeller {

    min_images = 10
      .type = int(value_min=1)
      .help = "The minimum number of images per pixel"

    filter_type = *median mean
      .type = choice
      .help = "The filter to use on the polar transformed image"

    kernel_size = 50
      .type = int(value_min=0)
      .help = "The kernel size for the median filter"

    niter = 100
      .type = int(value_min=1)
      .help = "The number of iterations for filling holes"

    image_type = min *mean
      .type = choice
      .help = "Which image to use"

  }

  include scope dials.algorithms.integration.integrator.phil_scope
  include scope dials.algorithms.spot_prediction.reflection_predictor.phil_scope

""",
    process_includes=True,
)


class ImageGenerator(object):
    """
    Generate diagnostic images

    """

    def __init__(self, model):
        """
        Init the model

        """
        import matplotlib

        matplotlib.use("Agg")
        self.model = model

    def save_min(self, filename):
        """
        Save the mean image

        """
        from matplotlib import pylab

        for i in range(len(self.model)):
            min_image = self.model[i].min_image
            vmax = sorted(list(min_image))[int(0.99 * len(min_image))]
            figure = pylab.figure(figsize=(6, 4))
            pylab.imshow(
                min_image.as_numpy_array(), interpolation="none", vmin=0, vmax=vmax
            )
            ax1 = pylab.gca()
            ax1.get_xaxis().set_visible(False)
            ax1.get_yaxis().set_visible(False)
            cb = pylab.colorbar()
            cb.ax.tick_params(labelsize=8)
            logger.info("Saving min image for panel %d to %s_%d.png" % (i, filename, i))
            pylab.savefig("%s_%d.png" % (filename, i), dpi=600, bbox_inches="tight")

    def save_max(self, filename):
        """
        Save the mean image

        """
        from matplotlib import pylab

        for i in range(len(self.model)):
            max_image = self.model[i].max_image
            vmax = sorted(list(max_image))[int(0.99 * len(max_image))]
            figure = pylab.figure(figsize=(6, 4))
            pylab.imshow(
                max_image.as_numpy_array(), interpolation="none", vmin=0, vmax=vmax
            )
            ax1 = pylab.gca()
            ax1.get_xaxis().set_visible(False)
            ax1.get_yaxis().set_visible(False)
            cb = pylab.colorbar()
            cb.ax.tick_params(labelsize=8)
            logger.info("Saving max image for panel %d to %s_%d.png" % (i, filename, i))
            pylab.savefig("%s_%d.png" % (filename, i), dpi=600, bbox_inches="tight")

    def save_mean(self, filename):
        """
        Save the mean image

        """
        from matplotlib import pylab

        for i in range(len(self.model)):
            mean = self.model[i].mean
            vmax = sorted(list(mean))[int(0.99 * len(mean))]
            figure = pylab.figure(figsize=(6, 4))
            pylab.imshow(mean.as_numpy_array(), interpolation="none", vmin=0, vmax=vmax)
            ax1 = pylab.gca()
            ax1.get_xaxis().set_visible(False)
            ax1.get_yaxis().set_visible(False)
            cb = pylab.colorbar()
            cb.ax.tick_params(labelsize=8)
            logger.info(
                "Saving mean image for panel %d to %s_%d.png" % (i, filename, i)
            )
            pylab.savefig("%s_%d.png" % (filename, i), dpi=600, bbox_inches="tight")

    def save_variance(self, filename):
        """
        Save the variance image

        """
        from matplotlib import pylab

        for i in range(len(self.model)):
            variance = self.model[i].variance
            vmax = sorted(list(variance))[int(0.99 * len(variance))]
            figure = pylab.figure(figsize=(6, 4))
            pylab.imshow(
                variance.as_numpy_array(), interpolation="none", vmin=0, vmax=vmax
            )
            ax1 = pylab.gca()
            ax1.get_xaxis().set_visible(False)
            ax1.get_yaxis().set_visible(False)
            cb = pylab.colorbar()
            cb.ax.tick_params(labelsize=8)
            logger.info(
                "Saving variance image for panel %d to %s_%d.png" % (i, filename, i)
            )
            pylab.savefig("%s_%d.png" % (filename, i), dpi=600, bbox_inches="tight")

    def save_dispersion(self, filename):
        """
        Save the dispersion image

        """
        from matplotlib import pylab

        for i in range(len(self.model)):
            dispersion = self.model[i].dispersion
            figure = pylab.figure(figsize=(6, 4))
            pylab.imshow(
                dispersion.as_numpy_array(), interpolation="none", vmin=0, vmax=2
            )
            ax1 = pylab.gca()
            ax1.get_xaxis().set_visible(False)
            ax1.get_yaxis().set_visible(False)
            cb = pylab.colorbar()
            cb.ax.tick_params(labelsize=8)
            logger.info(
                "Saving dispersion image for panel %d to %s_%d.png" % (i, filename, i)
            )
            pylab.savefig("%s_%d.png" % (filename, i), dpi=600, bbox_inches="tight")

    def save_mask(self, filename):
        """
        Save the dispersion image

        """
        from matplotlib import pylab

        for i in range(len(self.model)):
            mask = self.model[i].mask
            figure = pylab.figure(figsize=(6, 4))
            pylab.imshow(mask.as_numpy_array(), interpolation="none")
            ax1 = pylab.gca()
            ax1.get_xaxis().set_visible(False)
            ax1.get_yaxis().set_visible(False)
            logger.info(
                "Saving mask image for panel %d to %s_%d.png" % (i, filename, i)
            )
            pylab.savefig("%s_%d.png" % (filename, i), dpi=600, bbox_inches="tight")

    def save_model(self, filename):
        """
        Save the model image

        """
        from matplotlib import pylab

        for i in range(len(self.model)):
            model = self.model[i].model
            vmax = sorted(list(model))[int(0.99 * len(model))]
            figure = pylab.figure(figsize=(6, 4))
            pylab.imshow(
                model.as_numpy_array(), interpolation="none", vmin=0, vmax=vmax
            )
            ax1 = pylab.gca()
            ax1.get_xaxis().set_visible(False)
            ax1.get_yaxis().set_visible(False)
            cb = pylab.colorbar()
            cb.ax.tick_params(labelsize=8)
            logger.info(
                "Saving model image for panel %d to %s_%d.png" % (i, filename, i)
            )
            pylab.savefig("%s_%d.png" % (filename, i), dpi=600, bbox_inches="tight")

    def save_polar_model(self, filename):
        """
        Save the polar model image

        """
        from matplotlib import pylab

        for i in range(len(self.model)):
            polar_model = self.model[i].polar_model
            figure = pylab.figure(figsize=(6, 4))
            pylab.imshow(polar_model.as_numpy_array(), interpolation="none")
            ax1 = pylab.gca()
            ax1.get_xaxis().set_visible(False)
            ax1.get_yaxis().set_visible(False)
            cb = pylab.colorbar()
            cb.ax.tick_params(labelsize=8)
            logger.info(
                "Saving polar model image for panel %d to %s_%d.png" % (i, filename, i)
            )
            pylab.savefig("%s_%d.png" % (filename, i), dpi=600, bbox_inches="tight")


class Script(object):
    """A class for running the script."""

    def __init__(self):
        """Initialise the script."""
        from dials.util.options import OptionParser
        import libtbx.load_env

        # The script usage
        usage = (
            "usage: %s [options] [param.phil] "
            "models.expt" % libtbx.env.dispatcher_name
        )

        # Initialise the base class
        self.parser = OptionParser(
            usage=usage, phil=phil_scope, epilog=help_message, read_experiments=True
        )

    def run(self):
        """Execute the script."""
        from dials.util.command_line import heading
        from dials.array_family import flex
        from dials.util.options import flatten_experiments
        from time import time
        from dials.util import log
        from dials.util import Sorry
        from dials.algorithms.background.modeller import BackgroundModeller

        start_time = time()

        # Parse the command line
        params, options = self.parser.parse_args(show_diff_phil=False)

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

        # Only handle a single imageset at once
        imagesets = {expr.imageset for expr in experiments}
        if len(imagesets) != 1:
            raise Sorry("Can only process a single imageset at a time")

        # Predict the reflections
        logger.info("")
        logger.info("=" * 80)
        logger.info("")
        logger.info(heading("Predicting reflections"))
        logger.info("")
        predicted = flex.reflection_table.from_predictions_multi(
            experiments,
            dmin=params.prediction.d_min,
            dmax=params.prediction.d_max,
            margin=params.prediction.margin,
            force_static=params.prediction.force_static,
        )

        # Create the modeller
        modeller = BackgroundModeller(experiments, predicted, params)
        model = modeller.compute()

        # Save the background model
        logger.info("Saving background model to %s" % params.output.model)
        from dials.algorithms.background.gmodel import StaticBackgroundModel

        static_model = StaticBackgroundModel()
        for i in range(len(model)):
            static_model.add(model[i].model)
        with open(params.output.model, "wb") as outfile:
            import six.moves.cPickle as pickle

            pickle.dump(static_model, outfile, protocol=pickle.HIGHEST_PROTOCOL)

        # Output some diagnostic images
        image_generator = ImageGenerator(model)
        image_generator.save_mean(params.output.mean_image_prefix)
        image_generator.save_variance(params.output.variance_image_prefix)
        image_generator.save_dispersion(params.output.dispersion_image_prefix)
        image_generator.save_mask(params.output.mask_image_prefix)
        image_generator.save_min(params.output.min_image_prefix)
        image_generator.save_max(params.output.max_image_prefix)
        image_generator.save_model(params.output.model_image_prefix)
        # image_generator.save_polar_model(params.output.polar_model_image_prefix)

        # Print the time
        logger.info("Time Taken: %f" % (time() - start_time))


if __name__ == "__main__":
    from dials.util import halraiser

    try:
        script = Script()
        script.run()
    except Exception as e:
        halraiser(e)
