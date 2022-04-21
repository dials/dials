from __future__ import annotations

import logging
import pickle
import sys

from libtbx.phil import parse

import dials.util
import dials.util.log

logger = logging.getLogger("dials.command_line.model_background")

help_message = """
"""

# Set the phil scope
phil_scope = parse(
    """

  output {
    model = 'background.pickle'
      .type = str
      .help = "The output filename"

    log = 'dials.model_background.log'
      .type = str
      .help = "The log filename"

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


class ImageGenerator:
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

    def _save_plot(
        self, name, filename, extractor_fn, bounded=True, colorbar=True, vmax=None
    ):
        """
        Save the image
        """
        from matplotlib import pylab

        for i, model in enumerate(self.model):
            image = extractor_fn(model)
            pylab.figure(figsize=(6, 4))
            if bounded and vmax is None:
                boundaries = {
                    "vmin": 0,
                    "vmax": sorted(image)[int(0.99 * len(image))],
                }
            elif bounded:
                boundaries = {"vmin": 0, "vmax": vmax}
            else:
                boundaries = {}
            pylab.imshow(image.as_numpy_array(), interpolation="none", **boundaries)
            ax1 = pylab.gca()
            ax1.get_xaxis().set_visible(False)
            ax1.get_yaxis().set_visible(False)
            if colorbar:
                cb = pylab.colorbar()
                cb.ax.tick_params(labelsize=8)
            logger.info(
                "Saving %s image for panel %d to %s_%d.png", name, i, filename, i
            )
            pylab.savefig("%s_%d.png" % (filename, i), dpi=600, bbox_inches="tight")

    def save_min(self, filename):
        """
        Save the min image
        """
        self._save_plot("min", filename, lambda m: m.min_image)

    def save_max(self, filename):
        """
        Save the max image
        """
        self._save_plot("max", filename, lambda m: m.max_image)

    def save_mean(self, filename):
        """
        Save the mean image
        """
        self._save_plot("mean", filename, lambda m: m.mean)

    def save_variance(self, filename):
        """
        Save the variance image
        """
        self._save_plot("variance", filename, lambda m: m.variance)

    def save_dispersion(self, filename):
        """
        Save the dispersion image
        """
        self._save_plot("dispersion", filename, lambda m: m.dispersion, vmax=2)

    def save_mask(self, filename):
        """
        Save the dispersion image
        """
        self._save_plot(
            "mask", filename, lambda m: m.mask, bounded=False, colorbar=False
        )

    def save_model(self, filename):
        """
        Save the model image
        """
        self._save_plot("model", filename, lambda m: m.model)

    def save_polar_model(self, filename):
        """
        Save the polar model image
        """
        self._save_plot("polar model", filename, lambda m: m.polar_model, bounded=False)


class Script:
    """A class for running the script."""

    def __init__(self):
        """Initialise the script."""
        from dials.util.options import ArgumentParser

        usage = "dials.model_background [options] [param.phil] models.expt"

        # Initialise the base class
        self.parser = ArgumentParser(
            usage=usage, phil=phil_scope, epilog=help_message, read_experiments=True
        )

    def run(self, args=None):
        """Execute the script."""
        from dials.algorithms.background.modeller import BackgroundModeller
        from dials.array_family import flex
        from dials.util.command_line import heading
        from dials.util.options import flatten_experiments

        # Parse the command line
        params, options = self.parser.parse_args(args, show_diff_phil=False)

        # Configure the logging
        dials.util.log.config(verbosity=options.verbose, logfile=params.output.log)

        if params.integration.mp.nproc != 1 or params.integration.mp.njobs != 1:
            # https://github.com/dials/dials/issues/1083
            logger.warning(
                "Multiprocessing is currently disabled. " "Setting nproc = njobs = 1"
            )
            params.integration.mp.nproc = 1
            params.integration.mp.njobs = 1

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
            sys.exit("Can only process a single imageset at a time")

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
        logger.info("Saving background model to %s", params.output.model)
        from dials.algorithms.background.gmodel import StaticBackgroundModel

        static_model = StaticBackgroundModel()
        for m in model:
            static_model.add(m.model)
        with open(params.output.model, "wb") as outfile:
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


@dials.util.show_mail_handle_errors()
def run(args=None):
    script = Script()
    script.run(args)


if __name__ == "__main__":
    run()
