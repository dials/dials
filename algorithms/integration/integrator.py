from __future__ import absolute_import, division, print_function

import collections
import logging
import math
import random

import six.moves.cPickle as pickle
from dials_algorithms_integration_integrator_ext import *
from dials.algorithms.integration.processor import Processor3D
from dials.algorithms.integration.processor import ProcessorFlat3D
from dials.algorithms.integration.processor import Processor2D
from dials.algorithms.integration.processor import ProcessorSingle2D
from dials.algorithms.integration.processor import ProcessorStills
from dials.algorithms.integration.processor import ProcessorBuilder
from dials.algorithms.integration.processor import job
from dials.algorithms.integration.image_integrator import ImageIntegrator
from dials.array_family import flex
from dials.util import phil
from dials.util import Sorry

logger = logging.getLogger(__name__)


def generate_phil_scope():
    """
    Generate the integration phil scope.

    :return: The phil scope

    """
    import dials.extensions

    phil_scope = phil.parse(
        """

    integration {

      include scope dials.data.lookup.phil_scope

      block {

        size = auto
          .type = float
          .help = "The block size in rotation angle (degrees)."

        units = *degrees radians frames
          .type = choice
          .help = "The units of the block size"

        threshold = 0.99
          .type = float(value_min=0.0, value_max=1.0)
          .help = "For block size auto the block size is calculated by sorting"
                  "reflections by the number of frames they cover and then"
                  "selecting the block size to be 2*nframes[threshold] such"
                  "that 100*threshold % of reflections are guarenteed to be"
                  "fully contained in 1 block"

        force = False
          .type = bool
          .help = "If the number of processors is 1 and force is False, then the"
                  "number of blocks may be set to 1. If force is True then the"
                  "block size is always calculated."

        max_memory_usage = 0.75
          .type = float(value_min=0.0,value_max=1.0)
          .help = "The maximum percentage of total physical memory to use for"
                  "allocating shoebox arrays."

      }

      use_dynamic_mask = True
        .type = bool
        .help = "Use dynamic mask if available"

      debug {

        reference {

          filename = "reference_profiles.refl"
            .type = str
            .help = "The filename for the reference profiles"

          output = False
            .type = bool
            .help = "Save the reference profiles"
        }

        during = modelling *integration
          .type = choice
          .help = "Do debugging during modelling or integration"

        output = False
          .type = bool
          .help = "Save shoeboxes after each processing task."

        separate_files = True
          .type = bool
          .help = "If this is true, the shoeboxes are saved in separate files"
                  "from the output integrated.refl file. This is necessary"
                  "in most cases since the amount of memory used by the"
                  "shoeboxes is typically greater than the available system"
                  "memory. If, however, you know that memory is not an issue,"
                  "you can saved the shoeboxes in the integrated.refl file"
                  "by setting this option to False. This only works if the debug"
                  "output is during integrated and not modelling."

        delete_shoeboxes = False
          .type = bool
          .help = "Delete shoeboxes immediately before saving files. This option"
                  "in combination with debug.output=True enables intermediate"
                  "processing steps to make use of shoeboxes."

        select = None
          .type = reflection_table_selector
          .help = "A string specifying the selection. The string should be of the"
                  "form: select=${COLUMN}[<|<=|==|!=|>=|>]${VALUE}. In addition"
                  "to the items in the reflection table, the following implicit"
                  "columns are defined if the necessary data is there:"
                  " intensity.sum.i_over_sigma"
                  " intensity.prf.i_over_sigma"

        split_experiments = True
          .type = bool
          .help = "Split shoeboxes into different files"

      }

      integrator = *auto 3d flat3d 2d single2d stills volume 3d_threaded
        .type = choice
        .help = "The integrator to use."
        .expert_level=3

      profile {

        fitting = True
          .type = bool
          .help = "Use profile fitting if available"

        validation {

          number_of_partitions = 1
            .type = int(value_min=1)
            .help = "The number of subsamples to take from the reference spots."
                    "If the value is 1, then no validation is performed."

          min_partition_size = 100
            .type = int(value_min=1)
            .help = "The minimum number of spots to use in each subsample."

        }
      }

      filter
        .expert_level = 1
      {
        min_zeta = 0.05
          .help = "Filter the reflections by the value of zeta. A value of less"
                  "than or equal to zero indicates that this will not be used. A"
                  "positive value is used as the minimum permissable value."
          .type = float(value_min=0.0, value_max=1.0)

        max_shoebox_overlap = 1.0
          .type = float(value_min=0.0, value_max=1.0)
          .help = "Filter reflections whose shoeboxes are overlapped by greater"
                  "than the requested amount. Note that this is not the"
                  "percentage of the peak that is overlapped but rather the"
                  "percentage of the shoebox (background and foreground). This"
                  "can be useful when the detector is too close and many"
                  "overlapping reflections are predicted at high resolution"
                  "causing memory issues."

        ice_rings = False
          .help = "Set the ice ring flags"
          .type = bool
      }

      include scope dials.algorithms.integration.overlaps_filter.phil_scope

      mp {
        method = *none drmaa sge lsf pbs
          .type = choice
          .help = "The cluster method to use"

        njobs = 1
          .type = int(value_min=1)
          .help = "The number of cluster jobs to use"

        nproc = 1
          .type = int(value_min=1)
          .help = "The number of processes to use per cluster job"
      }

      summation {
        detector_gain = 1
          .type = float
          .help = "Multiplier for variances after integration of still images."
                  "See Leslie 1999."
      }
    }
  """,
        process_includes=True,
    )
    main_scope = phil_scope.get_without_substitution("integration")
    assert len(main_scope) == 1
    main_scope = main_scope[0]
    main_scope.adopt_scope(dials.extensions.Background.phil_scope())
    main_scope.adopt_scope(dials.extensions.Centroid.phil_scope())
    return phil_scope


# The integration phil scope
phil_scope = generate_phil_scope()


def hist(data, width=80, symbol="#", prefix=""):
    """
    A utility function to print a histogram of reflections on frames.

    :param data: The data to histogram
    :param width: The number of characters in each line
    :param symbol: The plot symbol
    :param prefix: String to prefix to each line
    :return: The histogram string

    """
    assert len(data) > 0, "Need > 0 reflections"
    assert width > 0, "Width should be > 0"
    count = collections.Counter(data)
    count = sorted(count.items())
    frame, count = zip(*count)
    max_frame = max(frame)
    min_count = min(count)
    max_count = max(count)
    assert max_count > 0, "Max should be > 0"
    assert min_count >= 0, "Min should be >= 0"
    if max_frame == 0:
        num_frame_zeros = 1
    else:
        num_frame_zeros = int(math.floor(math.log10(max_frame))) + 1
    num_count_zeros = int(math.floor(math.log10(max_count))) + 1
    assert num_frame_zeros > 0, "Num should be > 0"
    assert num_count_zeros > 0, "Num should be > 0"
    num_hist = width - (num_frame_zeros + num_count_zeros + 5) - len(prefix)
    assert num_hist > 0, "num_hist should be > 0"
    fmt = "%s%%-%dd [%%-%dd]: %%s" % (prefix, num_frame_zeros, num_count_zeros)
    scale = float(num_hist) / max_count
    return "\n".join(
        (
            fmt % (key, value, int(value * scale) * symbol)
            for key, value in zip(frame, count)
        )
    )


def frame_hist(bbox, width=80, symbol="#", prefix=""):
    """
    A utility function to print a histogram of reflections on frames.

    :param bbox: The bounding boxes
    :param width: The width of each line
    :param symbol: The histogram symbol
    :param prefix: A string to prefix to each line
    :return: The histogram string

    """
    return hist(
        [z for b in bbox for z in range(b[4], b[5])],
        width=width,
        symbol=symbol,
        prefix=prefix,
    )


def nframes_hist(bbox, width=80, symbol="#", prefix=""):
    """
    A utility function to print a histogram of number of frames.

    :param bbox: The bounding boxes
    :param width: The width of each line
    :param symbol: The histogram symbol
    :param prefix: A string to prefix to each line
    :return: The histogram string

    """
    return hist([b[5] - b[4] for b in bbox], width=width, symbol=symbol, prefix=prefix)


class Parameters(object):
    """
    A class to represent the integration parameters

    """

    class Filter(object):
        """
        Filter parameters

        """

        def __init__(self):
            self.min_zeta = 0.05
            self.powder_filter = None

    class Profile(object):
        """
        Profile parameters

        """

        class Validation(object):
            def __init__(self):
                self.number_of_partitions = 2
                self.min_partition_size = 100

        def __init__(self):
            self.fitting = True
            self.validation = Parameters.Profile.Validation()

    def __init__(self):
        """
        Initialize

        """
        from dials.algorithms.integration import processor

        self.modelling = processor.Parameters()
        self.integration = processor.Parameters()
        self.filter = Parameters.Filter()
        self.profile = Parameters.Profile()
        self.debug_reference_filename = "reference_profiles.refl"
        self.debug_reference_output = False

    @staticmethod
    def from_phil(params):
        """
        Convert the phil parameters

        """
        from dials.algorithms.integration import processor
        from dials.algorithms.integration.filtering import IceRingFilter

        # Init the parameters
        result = Parameters()

        # Create the multi processing parameters
        mp = processor.MultiProcessing()
        mp.method = params.mp.method
        mp.nproc = params.mp.nproc
        mp.njobs = params.mp.njobs

        # Set the lookup parameters
        lookup = processor.Lookup()
        lookup.mask = params.lookup.mask

        # Set the block parameters
        block = processor.Block()
        block.size = params.block.size
        block.units = params.block.units
        block.threshold = params.block.threshold
        block.force = params.block.force
        block.max_memory_usage = params.block.max_memory_usage

        # Set the modelling processor parameters
        result.modelling.mp = mp
        result.modelling.lookup = lookup
        result.modelling.block = block
        if params.debug.during == "modelling":
            result.modelling.debug.output = params.debug.output
        result.modelling.debug.select = params.debug.select
        result.modelling.debug.separate_files = True

        # Set the integration processor parameters
        result.integration.mp = mp
        result.integration.lookup = lookup
        result.integration.block = block
        if params.debug.during == "integration":
            result.integration.debug.output = params.debug.output
        result.integration.debug.select = params.debug.select
        result.integration.debug.separate_files = params.debug.separate_files
        result.integration.summation = params.summation

        result.debug_reference_filename = params.debug.reference.filename
        result.debug_reference_output = params.debug.reference.output

        # Get the min zeta filter
        result.filter.min_zeta = params.filter.min_zeta
        if params.filter.ice_rings is True:
            result.filter.powder_filter = IceRingFilter()

        # Get post-integration overlap filtering parameters
        result.integration.overlaps_filter = params.overlaps_filter

        # Set the profile fitting parameters
        result.profile.fitting = params.profile.fitting
        result.profile.validation.number_of_partitions = (
            params.profile.validation.number_of_partitions
        )
        result.profile.validation.min_partition_size = (
            params.profile.validation.min_partition_size
        )

        # Return the result
        return result


class InitializerRot(object):
    """
    A pre-processing class for oscillation data.

    """

    def __init__(self, experiments, params):
        """
        Initialise the pre-processor.

        """
        self.experiments = experiments
        self.params = params

    def __call__(self, reflections):
        """
        Do some pre-processing.
        """
        # Compute some reflection properties
        reflections.compute_zeta_multi(self.experiments)
        reflections.compute_d(self.experiments)
        reflections.compute_bbox(self.experiments)

        # Filter the reflections by zeta
        mask = flex.abs(reflections["zeta"]) < self.params.filter.min_zeta
        reflections.set_flags(mask, reflections.flags.dont_integrate)

        # Filter the reflections by powder ring
        if self.params.filter.powder_filter is not None:
            mask = self.params.filter.powder_filter(reflections["d"])
            reflections.set_flags(mask, reflections.flags.in_powder_ring)


class InitializerStills(object):
    """
    A pre-processing class for stills data.

    """

    def __init__(self, experiments, params):
        """
        Initialise the pre-processor.

        """
        self.experiments = experiments
        self.params = params

    def __call__(self, reflections):
        """
        Do some pre-processing.
        """
        # Compute some reflection properties
        reflections.compute_d(self.experiments)
        reflections.compute_bbox(self.experiments)

        # Check the bounding boxes are all 1 frame in width
        z0, z1 = reflections["bbox"].parts()[4:6]
        assert (z1 - z0).all_eq(1), "bbox is invalid"

        # Filter the reflections by powder ring
        if self.params.filter.powder_filter is not None:
            mask = self.params.filter.powder_filter(reflections["d"])
            reflections.set_flags(mask, reflections.flags.in_powder_ring)


class FinalizerBase(object):
    def __init__(self, reflections, experiments, params):
        """
        Initialise the post processor.

        """
        self.reflections = reflections
        self.experiments = experiments
        self.params = params

    def __call__(self):
        overlaps_scope = self.params.integration.overlaps_filter
        if True in [
            overlaps_scope.foreground_foreground.enable,
            overlaps_scope.foreground_background.enable,
        ]:
            from dials.algorithms.integration.overlaps_filter import (
                OverlapsFilterMultiExpt,
            )

            overlaps_filter = OverlapsFilterMultiExpt(
                self.reflections, self.experiments
            )
            if overlaps_scope.foreground_foreground.enable:
                overlaps_filter.remove_foreground_foreground_overlaps()
            if overlaps_scope.foreground_background.enable:
                overlaps_filter.remove_foreground_background_overlaps()
            self.reflections = overlaps_filter.refl


class FinalizerRot(FinalizerBase):
    """
    A post-processing class for oscillation data.

    """

    def __call__(self):
        """
        Do some post processing.

        """
        super(FinalizerRot, self).__call__()

        # Compute the corrections
        self.reflections.compute_corrections(self.experiments)


class FinalizerStills(FinalizerBase):
    """
    A post-processing class for stills data.

    """

    def __call__(self):
        """
        Do some post processing.

        """
        super(FinalizerStills, self).__call__()

        integrated = self.reflections

        # Select only those reflections which were integrated
        if "intensity.prf.variance" in integrated:
            selection = integrated.get_flags(integrated.flags.integrated, all=True)
        else:
            selection = integrated.get_flags(integrated.flags.integrated_sum)
        integrated = integrated.select(selection)

        len_all = len(integrated)
        integrated = integrated.select(
            ~integrated.get_flags(integrated.flags.foreground_includes_bad_pixels)
        )
        logger.info(
            "Filtering %d reflections with at least one bad foreground pixel out of %d"
            % (len_all - len(integrated), len_all)
        )

        # verify sigmas are sensible
        if "intensity.prf.value" in integrated:
            if (integrated["intensity.prf.variance"] <= 0).count(True) > 0:
                raise Sorry(
                    "Found negative variances (prf). Are bad pixels properly masked out?"
                )
        if "intensity.sum.value" in integrated:
            if (integrated["intensity.sum.variance"] <= 0).count(True) > 0:
                if (integrated["intensity.sum.variance"] < 0).count(True) > 0:
                    raise Sorry(
                        "Found negative variances (sum). Are bad pixels properly masked out?"
                    )
                n = (integrated["intensity.sum.variance"] == 0).count(True)
                sel = (integrated["intensity.sum.variance"] == 0) & (
                    integrated["intensity.sum.value"] == 0
                )
                if n == sel.count(True):
                    logger.info(
                        "Filtering %d reflections with no integrated signal (sum and variance = 0) out of %d"
                        % (n, len(integrated))
                    )
                    integrated = integrated.select(
                        integrated["intensity.sum.variance"] > 0
                    )
                else:
                    raise Sorry(
                        "Found reflections with variances == 0 but summed signal != 0"
                    )

            # apply detector gain to summation variances
            integrated[
                "intensity.sum.variance"
            ] *= self.params.integration.summation.detector_gain
        if "background.sum.value" in integrated:
            if (integrated["background.sum.variance"] < 0).count(True) > 0:
                raise Sorry(
                    "Found negative variances (background sum). Are bad pixels properly masked out?"
                )
            if (integrated["background.sum.variance"] == 0).count(True) > 0:
                logger.info(
                    "Filtering %d reflections with zero background variance"
                    % ((integrated["background.sum.variance"] == 0).count(True))
                )
                integrated = integrated.select(
                    integrated["background.sum.variance"] > 0
                )
            # apply detector gain to background summation variances
            integrated[
                "background.sum.variance"
            ] *= self.params.integration.summation.detector_gain

        self.reflections = integrated


class ProfileModellerExecutor(Executor):
    """
    The class to do profile modelling calculations

    """

    def __init__(self, experiments, profile_fitter):
        """
        Initialise the executor

        :param experiments: The experiment list

        """
        self.experiments = experiments
        self.profile_fitter = profile_fitter
        super(ProfileModellerExecutor, self).__init__()

    def initialize(self, frame0, frame1, reflections):
        """
        Initialise the processing for a job

        :param frame0: The first frame in the job
        :param frame1: The last frame in the job
        :param reflections: The reflections that will be processed

        """

        # Get some info
        EPS = 1e-7
        full_value = 0.997300203937 - EPS
        fully_recorded = reflections["partiality"] > full_value
        npart = fully_recorded.count(False)
        nfull = fully_recorded.count(True)
        nice = reflections.get_flags(reflections.flags.in_powder_ring).count(True)
        ntot = len(reflections)

        # Write some output
        logger.info(" Beginning modelling job %d" % job.index)
        logger.info("")
        logger.info(" Frames: %d -> %d" % (frame0, frame1))
        logger.info("")
        logger.info(" Number of reflections")
        logger.info("  Partial:     %d" % npart)
        logger.info("  Full:        %d" % nfull)
        logger.info("  In ice ring: %d" % nice)
        logger.info("  Total:       %d" % ntot)
        logger.info("")

        # Print a histogram of reflections on frames
        if frame1 - frame0 > 1:
            logger.info(
                " The following histogram shows the number of reflections predicted"
            )
            logger.info(" to have all or part of their intensity on each frame.")
            logger.info("")
            logger.info(frame_hist(reflections["bbox"], prefix=" ", symbol="*"))
            logger.info("")

    def process(self, frame, reflections):
        """
        Process the data

        :param frame: The frame being processed
        :param reflections: The reflections to process

        """
        # Check if pixels are overloaded
        reflections.is_overloaded(self.experiments)

        # Compute the shoebox mask
        reflections.compute_mask(self.experiments)

        # Process the data
        reflections.compute_background(self.experiments)
        reflections.compute_centroid(self.experiments)
        reflections.compute_summed_intensity()

        # Do the profile modelling
        self.profile_fitter.model(reflections)

        # Print some info
        fmt = " Modelled % 5d / % 5d reflection profiles on image %d"
        nmod = reflections.get_flags(reflections.flags.used_in_modelling).count(True)
        ntot = len(reflections)
        logger.info(fmt % (nmod, ntot, frame))

    def finalize(self):
        """
        Finalize the processing

        """
        pass

    def data(self):
        """
        :return: the modeller

        """
        return self.profile_fitter

    def __getinitargs__(self):
        """
        Support for pickling

        """
        return (self.experiments, self.profile_fitter)


class ProfileValidatorExecutor(Executor):
    """
    The class to do profile validation calculations

    """

    def __init__(self, experiments, profile_fitter):
        """
        Initialise the executor

        :param experiments: The experiment list

        """
        self.experiments = experiments
        self.profile_fitter = profile_fitter
        super(ProfileValidatorExecutor, self).__init__()

    def initialize(self, frame0, frame1, reflections):
        """
        Initialise the processing for a job

        :param frame0: The first frame in the job
        :param frame1: The last frame in the job
        :param reflections: The reflections that will be processed

        """

        # Get some info
        EPS = 1e-7
        full_value = 0.997300203937 - EPS
        fully_recorded = reflections["partiality"] > full_value
        npart = fully_recorded.count(False)
        nfull = fully_recorded.count(True)
        nice = reflections.get_flags(reflections.flags.in_powder_ring).count(True)
        ntot = len(reflections)

        # Write some output
        logger.info(" Beginning modelling job %d" % job.index)
        logger.info("")
        logger.info(" Frames: %d -> %d" % (frame0, frame1))
        logger.info("")
        logger.info(" Number of reflections")
        logger.info("  Partial:     %d" % npart)
        logger.info("  Full:        %d" % nfull)
        logger.info("  In ice ring: %d" % nice)
        logger.info("  Total:       %d" % ntot)
        logger.info("")

        # Print a histogram of reflections on frames
        if frame1 - frame0 > 1:
            logger.info(
                " The following histogram shows the number of reflections predicted"
            )
            logger.info(" to have all or part of their intensity on each frame.")
            logger.info("")
            logger.info(frame_hist(reflections["bbox"], prefix=" ", symbol="*"))
            logger.info("")

        self.results = None

    def process(self, frame, reflections):
        """
        Process the data

        :param frame: The frame being processed
        :param reflections: The reflections to process

        """
        # Check if pixels are overloaded
        reflections.is_overloaded(self.experiments)

        # Compute the shoebox mask
        reflections.compute_mask(self.experiments)

        # Process the data
        reflections.compute_background(self.experiments)
        reflections.compute_centroid(self.experiments)
        reflections.compute_summed_intensity()

        # Do the profile validation
        self.results = self.profile_fitter.validate(reflections)

        # Print some info
        fmt = " Validated % 5d / % 5d reflection profiles on image %d"
        nmod = reflections.get_flags(reflections.flags.used_in_modelling).count(True)
        ntot = len(reflections)
        logger.info(fmt % (nmod, ntot, frame))

    def finalize(self):
        """
        Finalize the processing

        """
        pass

    def data(self):
        """
        :return: the modeller

        """
        return self.results

    def __getinitargs__(self):
        """
        Support for pickling

        """
        return (self.experiments, self.profile_fitter)


class IntegratorExecutor(Executor):
    """
    The class to process the integration data

    """

    def __init__(self, experiments, profile_fitter=None):
        """
        Initialize the executor

        :param experiments: The experiment list

        """
        self.experiments = experiments
        self.overlaps = None
        self.profile_fitter = profile_fitter
        super(IntegratorExecutor, self).__init__()

    def initialize(self, frame0, frame1, reflections):
        """
        Initialize the processing for the job

        :param frame0: The first frame to process
        :param frame1: The last frame to process
        :param reflections: The reflections to process

        """

        # Get some info
        EPS = 1e-7
        full_value = 0.997300203937 - EPS
        fully_recorded = reflections["partiality"] > full_value
        npart = fully_recorded.count(False)
        nfull = fully_recorded.count(True)
        nice = reflections.get_flags(reflections.flags.in_powder_ring).count(True)
        nint = reflections.get_flags(reflections.flags.dont_integrate).count(False)
        ntot = len(reflections)

        # Write some output
        logger.info(" Beginning integration job %d" % job.index)
        logger.info("")
        logger.info(" Frames: %d -> %d" % (frame0, frame1))
        logger.info("")
        logger.info(" Number of reflections")
        logger.info("  Partial:     %d" % npart)
        logger.info("  Full:        %d" % nfull)
        logger.info("  In ice ring: %d" % nice)
        logger.info("  Integrate:   %d" % nint)
        logger.info("  Total:       %d" % ntot)
        logger.info("")

        # Print a histogram of reflections on frames
        if frame1 - frame0 > 1:
            logger.info(
                " The following histogram shows the number of reflections predicted"
            )
            logger.info(" to have all or part of their intensity on each frame.")
            logger.info("")
            logger.info(frame_hist(reflections["bbox"], prefix=" ", symbol="*"))
            logger.info("")

        # Find any overlaps
        self.overlaps = reflections.find_overlaps(self.experiments)

    def process(self, frame, reflections):
        """
        Process the reflections on a frame

        :param frame: The frame to process
        :param reflections: The reflections to process

        """
        from dials.algorithms.shoebox import MaskCode

        # Check if pixels are overloaded
        reflections.is_overloaded(self.experiments)

        # Compute the shoebox mask
        reflections.compute_mask(self.experiments)

        # Check for invalid pixels in foreground/background
        reflections.contains_invalid_pixels()

        # Process the data
        reflections.compute_background(self.experiments)
        reflections.compute_centroid(self.experiments)
        reflections.compute_summed_intensity()
        if self.profile_fitter:
            reflections.compute_fitted_intensity(self.profile_fitter)

        # from matplotlib import pylab
        # from dials.array_family import flex
        # I = reflections['intensity.sum.value']
        # F = reflections.get_flags(reflections.flags.integrated_sum)
        # A = flex.size_t(range(len(reflections)))
        # A = A.select(F)
        # I = I.select(F)
        # index = flex.max_index(I)
        # index = A[index]
        # sbox = reflections['shoebox'][index]
        # zmin, zmax = sbox.bbox[4:6]
        # zind = int((zmin + zmax) / 2)
        # image1 = sbox.background.as_numpy_array().sum(axis=0)
        # image2 = sbox.data.as_numpy_array().sum(axis=0)
        # vmin=image1.min()#min([image1.min(), image2.min()])
        # vmax=image1.max()#min([image1.max(), image2.max()])
        # print image1.min(), image1.max(), image2.min(), image2.max()
        # pylab.subplot(121)
        # pylab.imshow(image1, interpolation='none', vmin=vmin, vmax=vmax)
        # pylab.subplot(122)
        # pylab.imshow(image2, interpolation='none', vmin=vmin, vmax=vmax)
        # pylab.show()

        # Compute the number of background/foreground pixels
        sbox = reflections["shoebox"]
        code1 = MaskCode.Valid
        code2 = MaskCode.Background | code1
        code3 = MaskCode.BackgroundUsed | code2
        code4 = MaskCode.Foreground | code1
        reflections["num_pixels.valid"] = sbox.count_mask_values(code1)
        reflections["num_pixels.background"] = sbox.count_mask_values(code2)
        reflections["num_pixels.background_used"] = sbox.count_mask_values(code3)
        reflections["num_pixels.foreground"] = sbox.count_mask_values(code4)

        # Print some info
        fmt = " Integrated % 5d (sum) + % 5d (prf) /% 5d reflections on image %d"
        nsum = reflections.get_flags(reflections.flags.integrated_sum).count(True)
        nprf = reflections.get_flags(reflections.flags.integrated_prf).count(True)
        ntot = len(reflections)
        logger.info(fmt % (nsum, nprf, ntot, frame))

    def finalize(self):
        """
        Finalize the processing

        """
        pass

    def data(self):
        """
        Return data

        """
        return None

    def __getinitargs__(self):
        """
        Support for pickling

        """
        return (self.experiments, self.profile_fitter)


class Integrator(object):
    """
    The integrator class

    """

    def __init__(self, experiments, reflections, params):
        """
        Initialize the integrator

        :param experiments: The experiment list
        :param reflections: The reflections to process
        :param params: The parameters to use

        """

        # Save some stuff
        self.experiments = experiments
        self.reflections = reflections
        self.params = Parameters.from_phil(params.integration)
        self.profile_model_report = None
        self.integration_report = None

    def integrate(self):
        """
        Integrate the data

        """
        from dials.algorithms.integration.report import IntegrationReport
        from dials.algorithms.integration.report import ProfileModelReport
        from dials.algorithms.integration.report import ProfileValidationReport
        from dials.util.command_line import heading
        from dials.util import pprint
        from dials.algorithms.profile_model.modeller import MultiExpProfileModeller
        from dials.algorithms.integration.validation import (
            ValidatedMultiExpProfileModeller,
        )

        # Ensure we get the same random sample each time
        random.seed(0)

        # Init the report
        self.profile_model_report = None
        self.integration_report = None

        # Heading
        logger.info("=" * 80)
        logger.info("")
        logger.info(heading("Processing reflections"))
        logger.info("")

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

        # Initialize the reflections
        initialize = self.InitializerClass(self.experiments, self.params)
        initialize(self.reflections)

        # Check if we want to do some profile fitting
        fitting_class = [e.profile.fitting_class() for e in self.experiments]
        fitting_avail = all(c is not None for c in fitting_class)
        if self.params.profile.fitting and fitting_avail:
            profile_fitting = True
            profile_fitter = None
        else:
            profile_fitting = False
            profile_fitter = None

        # Do profile modelling
        if profile_fitting:

            logger.info("=" * 80)
            logger.info("")
            logger.info(heading("Modelling reflection profiles"))
            logger.info("")

            # Get the selection
            selection = self.reflections.get_flags(
                self.reflections.flags.reference_spot
            )

            # Get the reference spots
            reference = self.reflections.select(selection)

            # Check if we need to skip
            if len(reference) == 0:
                logger.info(
                    "** Skipping profile modelling - no reference profiles given **"
                )
            else:

                # Try to set up the validation
                if self.params.profile.validation.number_of_partitions > 1:
                    n = len(reference)
                    k_max = int(
                        math.floor(
                            n / self.params.profile.validation.min_partition_size
                        )
                    )
                    if k_max < self.params.profile.validation.number_of_partitions:
                        num_folds = k_max
                    else:
                        num_folds = self.params.profile.validation.number_of_partitions
                    if num_folds > 1:
                        indices = (
                            list(range(num_folds)) * int(math.ceil(n / num_folds))
                        )[0:n]
                        random.shuffle(indices)
                        reference["profile.index"] = flex.size_t(indices)
                    if num_folds < 1:
                        num_folds = 1
                else:
                    num_folds = 1

                # Create the profile fitter
                profile_fitter = ValidatedMultiExpProfileModeller()
                for i in range(num_folds):
                    profile_fitter_single = MultiExpProfileModeller()  # (num_folds)
                    for expr in self.experiments:
                        profile_fitter_single.add(expr.profile.fitting_class()(expr))
                    profile_fitter.add(profile_fitter_single)

                # Create the data processor
                executor = ProfileModellerExecutor(self.experiments, profile_fitter)
                processor = ProcessorBuilder(
                    self.ProcessorClass,
                    self.experiments,
                    reference,
                    self.params.modelling,
                ).build()
                processor.executor = executor

                # Process the reference profiles
                reference, profile_fitter_list, time_info = processor.process()

                # Set the reference spots info
                # self.reflections.set_selected(selection, reference)

                # Finalize the profile models for validation
                assert len(profile_fitter_list) > 0, "No profile fitters"
                profile_fitter = None
                for pf in profile_fitter_list.values():
                    if pf is None:
                        continue
                    if profile_fitter is None:
                        profile_fitter = pf
                    else:
                        profile_fitter.accumulate(pf)
                profile_fitter.finalize()

                # Get the finalized modeller
                finalized_profile_fitter = profile_fitter.finalized_model()

                # Print profiles
                if self.params.debug_reference_output:
                    reference_debug = []
                    for i in range(len(finalized_profile_fitter)):
                        m = finalized_profile_fitter[i]
                        p = []
                        for j in range(len(m)):
                            try:
                                p.append((m.data(j), m.mask(j)))
                            except Exception:
                                p.append(None)
                    reference_debug.append(p)
                    with open(self.params.debug_reference_filename, "wb") as outfile:
                        pickle.dump(reference_debug, outfile)

                for i in range(len(finalized_profile_fitter)):
                    m = finalized_profile_fitter[i]
                    logger.debug("")
                    logger.debug("Profiles for experiment %d" % i)
                    for j in range(len(m)):
                        logger.debug("Profile %d" % j)
                        try:
                            logger.debug(pprint.profile3d(m.data(j)))
                        except Exception:
                            logger.debug("** NO PROFILE **")

                # Print the modeller report
                self.profile_model_report = ProfileModelReport(
                    self.experiments, finalized_profile_fitter, reference
                )
                logger.info("")
                logger.info(self.profile_model_report.as_str(prefix=" "))

                # Print the time info
                logger.info("")
                logger.info("Timing information for reference profile formation")
                logger.info(str(time_info))
                logger.info("")

                # If we have more than 1 fold then do the validation
                if num_folds > 1:

                    # Create the data processor
                    executor = ProfileValidatorExecutor(
                        self.experiments, profile_fitter
                    )
                    processor = ProcessorBuilder(
                        self.ProcessorClass,
                        self.experiments,
                        reference,
                        self.params.modelling,
                    ).build()
                    processor.executor = executor

                    # Process the reference profiles
                    reference, validation, time_info = processor.process()

                    # Print the modeller report
                    self.profile_validation_report = ProfileValidationReport(
                        self.experiments, profile_fitter, reference, num_folds
                    )
                    logger.info("")
                    logger.info(self.profile_validation_report.as_str(prefix=" "))

                    # Print the time info
                    logger.info("")
                    logger.info("Timing information for reference profile validation")
                    logger.info(str(time_info))
                    logger.info("")

                # Set to the finalized fitter
                profile_fitter = finalized_profile_fitter

        logger.info("=" * 80)
        logger.info("")
        logger.info(heading("Integrating reflections"))
        logger.info("")

        # Create the data processor
        executor = IntegratorExecutor(self.experiments, profile_fitter)
        processor = ProcessorBuilder(
            self.ProcessorClass,
            self.experiments,
            self.reflections,
            self.params.integration,
        ).build()
        processor.executor = executor

        # Process the reflections
        self.reflections, _, time_info = processor.process()

        # Finalize the reflections
        finalize = self.FinalizerClass(self.reflections, self.experiments, self.params)
        finalize()
        self.reflections = finalize.reflections
        self.experiments = finalize.experiments

        # Create the integration report
        self.integration_report = IntegrationReport(self.experiments, self.reflections)
        logger.info("")
        logger.info(self.integration_report.as_str(prefix=" "))

        # Print the time info
        logger.info("Timing information for integration")
        logger.info(str(time_info))
        logger.info("")

        # Return the reflections
        return self.reflections

    def report(self):
        """
        Return the report of the processing

        """
        from dials.util.report import Report

        result = Report()
        if self.profile_model_report is not None:
            result.combine(self.profile_model_report)
        result.combine(self.integration_report)
        return result

    def summary(self, block_size, block_size_units):
        """ Print a summary of the integration stuff. """
        from libtbx.table_utils import format as table

        # Compute the task table
        if self._experiments.all_stills():
            rows = [["#", "Group", "Frame From", "Frame To"]]
            for i in range(len(self)):
                job = self._manager.job(i)
                group = job.index()
                f0, f1 = job.frames()
                rows.append([str(i), str(group), str(f0), str(f1)])
        elif self._experiments.all_sweeps():
            rows = [["#", "Group", "Frame From", "Frame To", "Angle From", "Angle To"]]
            for i in range(len(self)):
                job = self._manager.job(i)
                group = job.index()
                expr = job.expr()
                f0, f1 = job.frames()
                scan = self._experiments[expr[0]].scan
                p0 = scan.get_angle_from_array_index(f0)
                p1 = scan.get_angle_from_array_index(f1)
                rows.append([str(i), str(group), str(f0), str(f1), str(p0), str(p1)])
        else:
            raise RuntimeError("Experiments must be all sweeps or all stills")
        task_table = table(rows, has_header=True, justify="right", prefix=" ")


class Integrator3D(Integrator):
    """
    Integrator for 3D algorithms

    """

    InitializerClass = InitializerRot
    ProcessorClass = Processor3D
    FinalizerClass = FinalizerRot


class IntegratorFlat3D(Integrator):
    """
    Integrator for flattened 3D algorithms

    """

    InitializerClass = InitializerRot
    ProcessorClass = ProcessorFlat3D
    FinalizerClass = FinalizerRot


class Integrator2D(Integrator):
    """
    Integrator for 2D algorithms

    """

    InitializerClass = InitializerRot
    ProcessorClass = Processor2D
    FinalizerClass = FinalizerRot


class IntegratorSingle2D(Integrator):
    """
    Integrator for 2D algorithms on a single image

    """

    InitializerClass = InitializerRot
    ProcessorClass = ProcessorSingle2D
    FinalizerClass = FinalizerRot


class IntegratorStills(Integrator):
    """
    Integrator for still algorithms

    """

    InitializerClass = InitializerStills
    ProcessorClass = ProcessorStills
    FinalizerClass = FinalizerStills


class IntegratorVolume(ImageIntegrator):
    """
    Volume integrator

    """

    pass


class Integrator3DThreaded(object):
    """
    Integrator for 3D algorithms

    """

    def __init__(self, experiments, reflections, params):

        """
        Initialize the integrator

        :param experiments: The experiment list
        :param reflections: The reflections to process
        :param params: The parameters to use

        """

        # Save some stuff
        self.experiments = experiments
        self.reflections = reflections
        self.params = params
        self.profile_model_report = None
        self.integration_report = None

    def initialise(self):
        """
        Initialise the integrator

        """
        # Compute some reflection properties
        self.reflections.compute_zeta_multi(self.experiments)
        self.reflections.compute_d(self.experiments)
        self.reflections.compute_bbox(self.experiments)

        # Filter the reflections by zeta
        mask = (
            flex.abs(self.reflections["zeta"]) < self.params.integration.filter.min_zeta
        )
        self.reflections.set_flags(mask, self.reflections.flags.dont_integrate)

        # Filter the reflections by powder ring
        # if self.params.integration.filter.powder_filter is not None:
        #   mask = self.params.integration.filter.powder_filter(self.reflections['d'])
        #   self.reflections.set_flags(mask, self.reflections.flags.in_powder_ring)

    def finalise(self):
        """
        Finalise the integrator

        """

        # Compute the corrections
        self.reflections.compute_corrections(self.experiments)

    def integrate(self):
        """
        Integrate the data

        """
        from dials.algorithms.integration.parallel_integrator import (
            ReferenceCalculatorProcessor,
        )
        from dials.algorithms.integration.parallel_integrator import IntegratorProcessor
        from dials.algorithms.integration.report import IntegrationReport
        from dials.util.command_line import heading

        # Init the report
        self.profile_model_report = None
        self.integration_report = None

        # Heading
        logger.info("=" * 80)
        logger.info("")
        logger.info(heading("Processing reflections"))
        logger.info("")

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

        # Do the initialisation
        self.initialise()

        # Do profile modelling
        if self.params.integration.profile.fitting:

            logger.info("=" * 80)
            logger.info("")
            logger.info(heading("Modelling reflection profiles"))
            logger.info("")

            # Compute the reference profiles
            reference_calculator = ReferenceCalculatorProcessor(
                experiments=self.experiments,
                reflections=self.reflections,
                params=self.params,
            )

            # Get the reference profiles
            self.reference_profiles = reference_calculator.profiles()
        else:
            self.reference_profiles = None

        logger.info("=" * 80)
        logger.info("")
        logger.info(heading("Integrating reflections"))
        logger.info("")

        integrator = IntegratorProcessor(
            experiments=self.experiments,
            reflections=self.reflections,
            reference=self.reference_profiles,
            params=self.params,
        )

        # Process the reflections
        self.reflections = integrator.reflections()

        # Do the finalisation
        self.finalise()

        # Create the integration report
        self.integration_report = IntegrationReport(self.experiments, self.reflections)
        logger.info("")
        logger.info(self.integration_report.as_str(prefix=" "))

        # Print the time info
        # logger.info("Timing information for integration")
        # logger.info(str(time_info))
        # logger.info("")

        # Return the reflections
        return self.reflections

    def report(self):
        """
        Return the report of the processing

        """
        from dials.util.report import Report

        result = Report()
        if self.profile_model_report is not None:
            result.combine(self.profile_model_report)
        result.combine(self.integration_report)
        return result

    def summary(self, block_size, block_size_units):
        """ Print a summary of the integration stuff. """
        from libtbx.table_utils import format as table

        # Compute the task table
        if self._experiments.all_stills():
            rows = [["#", "Group", "Frame From", "Frame To"]]
            for i in range(len(self)):
                job = self._manager.job(i)
                group = job.index()
                f0, f1 = job.frames()
                rows.append([str(i), str(group), str(f0), str(f1)])
        elif self._experiments.all_sweeps():
            rows = [["#", "Group", "Frame From", "Frame To", "Angle From", "Angle To"]]
            for i in range(len(self)):
                job = self._manager.job(i)
                group = job.index()
                expr = job.expr()
                f0, f1 = job.frames()
                scan = self._experiments[expr[0]].scan
                p0 = scan.get_angle_from_array_index(f0)
                p1 = scan.get_angle_from_array_index(f1)
                rows.append([str(i), str(group), str(f0), str(f1), str(p0), str(p1)])
        else:
            raise RuntimeError("Experiments must be all sweeps or all stills")
        task_table = table(rows, has_header=True, justify="right", prefix=" ")


class IntegratorFactory(object):
    """
    A factory for creating integrators.
    """

    @staticmethod
    def create(params, experiments, reflections):
        """
        Create the integrator from the input configuration.

        :param params: The input phil parameters
        :param experiments: The list of experiments
        :param reflections: The reflections to integrate
        :return: The integrator class
        """
        import dials.extensions
        from dials.util import Sorry

        # Check each experiment has an imageset
        for exp in experiments:
            if exp.imageset is None:
                raise Sorry(
                    """
          One or more experiment does not contain an imageset. Access to the
          image data is crucial for integration.
        """
                )

        # Read the mask in if necessary
        if params.integration.lookup.mask is not None:
            if isinstance(params.integration.lookup.mask, str):
                with open(params.integration.lookup.mask, "rb") as infile:
                    params.integration.lookup.mask = pickle.load(infile)

        # Initialise the strategy classes
        BackgroundAlgorithm = dials.extensions.Background.load(
            params.integration.background.algorithm
        )
        CentroidAlgorithm = dials.extensions.Centroid.load(
            params.integration.centroid.algorithm
        )

        # Set the algorithms in the reflection table
        flex.reflection_table._background_algorithm = flex.strategy(
            BackgroundAlgorithm, params
        )
        flex.reflection_table._centroid_algorithm = flex.strategy(
            CentroidAlgorithm, params
        )

        # Get the classes we need
        if params.integration.integrator == "auto":
            if experiments.all_stills():
                params.integration.integrator = "stills"
            else:
                params.integration.integrator = "3d"
        if params.integration.integrator == "3d":
            IntegratorClass = Integrator3D
        elif params.integration.integrator == "flat3d":
            IntegratorClass = IntegratorFlat3D
        elif params.integration.integrator == "2d":
            IntegratorClass = Integrator2D
        elif params.integration.integrator == "single2d":
            IntegratorClass = IntegratorSingle2D
        elif params.integration.integrator == "stills":
            IntegratorClass = IntegratorStills
        elif params.integration.integrator == "volume":
            IntegratorClass = IntegratorVolume
        elif params.integration.integrator == "3d_threaded":
            IntegratorClass = Integrator3DThreaded
        else:
            raise RuntimeError("Unknown integration type")

        # Remove scan if stills
        if experiments.all_stills():
            for experiment in experiments:
                experiment.scan = None

        # Return an instantiation of the class
        return IntegratorClass(experiments, reflections, params)
