from __future__ import absolute_import, division, print_function

import logging

from dials_algorithms_background_modeller_ext import *

logger = logging.getLogger(__name__)


class FinalizeModel(object):
    """
    A class to finalize the background model

    """

    def __init__(self, experiments, filter_type="median", kernel_size=10, niter=100):
        """
        Initialize the finalizer

        :param experiments: The experiment list
        :param kernel_size: The median filter kernel size
        :param niter: The number of iterations for filling holes

        """
        from dials.algorithms.background.gmodel import PolarTransform

        # Set some parameters
        self.filter_type = filter_type
        self.kernel_size = kernel_size
        self.niter = niter

        # Check the input
        assert len(experiments) == 1
        experiment = experiments[0]
        assert len(experiment.detector) == 1

        # Save the experiment
        self.experiment = experiment

        # Create the transform object
        self.transform = PolarTransform(
            experiment.beam, experiment.detector[0], experiment.goniometer
        )

    def finalize(self, data, mask):
        """
        Finalize the model

        :param data: The data array
        :param mask: The mask array

        """
        from dials.algorithms.image.filter import median_filter, mean_filter
        from dials.algorithms.image.fill_holes import diffusion_fill
        from dials.algorithms.image.fill_holes import simple_fill
        from dials.array_family import flex

        # Print some image properties
        sub_data = data.as_1d().select(mask.as_1d())
        logger.info("Raw image statistics:")
        logger.info("  min:  %d" % int(flex.min(sub_data)))
        logger.info("  max:  %d" % int(flex.max(sub_data)))
        logger.info("  mean: %d" % int(flex.mean(sub_data)))
        logger.info("")

        # Transform to polar
        logger.info("Transforming image data to polar grid")
        result = self.transform.to_polar(data, mask)
        data = result.data()
        mask = result.mask()
        sub_data = data.as_1d().select(mask.as_1d())
        logger.info("Polar image statistics:")
        logger.info("  min:  %d" % int(flex.min(sub_data)))
        logger.info("  max:  %d" % int(flex.max(sub_data)))
        logger.info("  mean: %d" % int(flex.mean(sub_data)))
        logger.info("")

        # Filter the image to remove noise
        if self.kernel_size > 0:
            if self.filter_type == "median":
                logger.info("Applying median filter")
                data = median_filter(data, mask, (self.kernel_size, 0), periodic=True)
                sub_data = data.as_1d().select(mask.as_1d())
                logger.info("Median polar image statistics:")
                logger.info("  min:  %d" % int(flex.min(sub_data)))
                logger.info("  max:  %d" % int(flex.max(sub_data)))
                logger.info("  mean: %d" % int(flex.mean(sub_data)))
                logger.info("")
            elif self.filter_type == "mean":
                logger.info("Applying mean filter")
                mask_as_int = mask.as_1d().as_int()
                mask_as_int.reshape(mask.accessor())
                data = mean_filter(data, mask_as_int, (self.kernel_size, 0), 1)
                sub_data = data.as_1d().select(mask.as_1d())
                logger.info("Mean polar image statistics:")
                logger.info("  min:  %d" % int(flex.min(sub_data)))
                logger.info("  max:  %d" % int(flex.max(sub_data)))
                logger.info("  mean: %d" % int(flex.mean(sub_data)))
                logger.info("")
            else:
                raise RuntimeError("Unknown filter_type: %s" % self.filter_type)

        # Fill any remaining holes
        logger.info("Filling holes")
        data = simple_fill(data, mask)
        data = diffusion_fill(data, mask, self.niter)
        mask = flex.bool(data.accessor(), True)
        sub_data = data.as_1d().select(mask.as_1d())
        logger.info("Filled polar image statistics:")
        logger.info("  min:  %d" % int(flex.min(sub_data)))
        logger.info("  max:  %d" % int(flex.max(sub_data)))
        logger.info("  mean: %d" % int(flex.mean(sub_data)))
        logger.info("")

        # Transform back
        logger.info("Transforming image data from polar grid")
        result = self.transform.from_polar(data, mask)
        data = result.data()
        mask = result.mask()
        sub_data = data.as_1d().select(mask.as_1d())
        logger.info("Final image statistics:")
        logger.info("  min:  %d" % int(flex.min(sub_data)))
        logger.info("  max:  %d" % int(flex.max(sub_data)))
        logger.info("  mean: %d" % int(flex.mean(sub_data)))
        logger.info("")

        # Fill in any discontinuities
        mask = ~self.transform.discontinuity()[:-1, :-1]
        data = diffusion_fill(data, mask, self.niter)

        # Get and apply the mask
        mask = self.experiment.imageset.get_mask(0)[0]
        mask = mask.as_1d().as_int().as_double()
        mask.reshape(data.accessor())
        data *= mask

        # Return the result
        return data


class BackgroundModellerResult(object):
    """
    A class to contain the modelling result

    """

    def __init__(
        self,
        mean=None,
        variance=None,
        dispersion=None,
        mask=None,
        min_image=None,
        max_image=None,
        model=None,
        polar_model=None,
    ):
        """
        Init the result

        """
        self.mean = mean
        self.variance = variance
        self.dispersion = dispersion
        self.mask = mask
        self.min_image = min_image
        self.max_image = max_image
        self.model = model
        self.polar_model = polar_model


class BackgroundModellerExecutor(object):
    def __init__(self, experiments, params):
        assert len(experiments) == 1
        self.min_images = params.modeller.min_images
        if self.min_images > len(experiments[0].imageset):
            self.min_images = len(experiments[0].imageset)
        self.image_type = params.modeller.image_type

        self.finalizer = FinalizeModel(
            experiments=experiments,
            filter_type=params.modeller.filter_type,
            kernel_size=params.modeller.kernel_size,
            niter=params.modeller.niter,
        )
        self.result = None

    def process(self, image_volume, experiments, reflections):
        from dials.algorithms.integration.processor import job

        # Write some output
        logger.info(
            " Background modelling; job: %d; frames: %d -> %d; # Reflections: %d"
            % (
                job.index,
                image_volume.frame0(),
                image_volume.frame1(),
                len(reflections),
            )
        )

        # Compute the shoebox mask
        reflections.compute_mask(experiments=experiments, image_volume=image_volume)

        # Compute the sum, sum^2 and the number of contributing pixels
        return MultiPanelBackgroundStatistics(image_volume)

    def accumulate(self, index, data):
        if self.result is None:
            self.result = data
        else:
            self.result += data

    def finalize_model(self):
        logger.info("")
        logger.info("=" * 80)
        logger.info("Finalizing model")
        logger.info("")

        result = []
        for i in range(len(self.result)):

            # Get the statistics
            stats = self.result.get(i)
            mean = stats.mean(self.min_images)
            variance = stats.variance(self.min_images)
            dispersion = stats.dispersion(self.min_images)
            mask = stats.mask(self.min_images)
            min_image = stats.min()
            max_image = stats.max()

            # Create the model
            if self.image_type == "min":
                model = self.finalizer.finalize(min_image, mask)
            elif self.image_type == "mean":
                model = self.finalizer.finalize(mean, mask)
            else:
                raise RuntimeError("Unknown image_type: %s" % self.image_type)

            # Add to the list
            result.append(
                BackgroundModellerResult(
                    mean=mean,
                    variance=variance,
                    dispersion=dispersion,
                    mask=mask,
                    min_image=min_image,
                    max_image=max_image,
                    model=model,
                )
            )

        return result


class BackgroundModeller(object):
    """
    A class to help with background modelling

    """

    def __init__(self, experiments, reflections, params):
        """
        Initialize the modeller

        :param experiments: The experiment list
        :param reflections: The reflections to process
        :param params: The parameters to use

        """
        # Check all reflections have same imageset and get it
        imageset = experiments[0].imageset
        for expr in experiments:
            assert expr.imageset == imageset, "All experiments must share and imageset"

        # Save some stuff
        self.experiments = experiments
        self.reflections = reflections
        self.params = params
        self.model = None

    def compute(self):
        """
        Integrate the data

        """
        from dials.algorithms.integration.image_integrator import ProcessorImage
        from dials.util.command_line import heading

        # Init the report
        self.profile_model_report = None
        self.integration_report = None

        # Create summary format
        fmt = (
            " Processing the following experiments:\n"
            "\n"
            " Experiments: %d\n"
            " Beams:       %d\n"
            " Detectors:   %d\n"
            " Goniometers: %d\n"
            " Scans:       %d\n"
            " Crystals:    %d\n"
            " Imagesets:   %d\n"
        )

        # Print the summary
        logger.info(
            fmt
            % (
                len(self.experiments),
                len(self.experiments.beams()),
                len(self.experiments.detectors()),
                len(self.experiments.goniometers()),
                len(self.experiments.scans()),
                len(self.experiments.crystals()),
                len(self.experiments.imagesets()),
            )
        )

        # Print a heading
        logger.info("=" * 80)
        logger.info("")
        logger.info(heading("Modelling background"))
        logger.info("")

        # Expand n_sigma
        for expt in self.experiments:
            expt.profile._n_sigma += 2

        # Compute some reflection properties
        self.reflections.compute_zeta_multi(self.experiments)
        self.reflections.compute_d(self.experiments)
        self.reflections.compute_bbox(self.experiments)

        # Construvt the image integrator processor
        processor = ProcessorImage(self.experiments, self.reflections, self.params)
        processor.executor = BackgroundModellerExecutor(self.experiments, self.params)

        # Do the processing
        _, time_info = processor.process()

        # Compute the model
        self.model = processor.executor.finalize_model()

        # Print the time info
        logger.info(str(time_info))
        logger.info("")

        # Return the reflections
        return self.model
