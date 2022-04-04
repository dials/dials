from __future__ import annotations

import collections
import functools
import logging
import math
import pickle
import random

import dials.extensions
from dials.algorithms.integration import TimingInfo, processor
from dials.algorithms.integration.filtering import IceRingFilter
from dials.algorithms.integration.parallel_integrator import (
    IntegratorProcessor,
    ReferenceCalculatorProcessor,
)
from dials.algorithms.integration.processor import (
    Processor2D,
    Processor3D,
    ProcessorFlat3D,
    ProcessorSingle2D,
    ProcessorStills,
    assess_available_memory,
    build_processor,
    job,
)
from dials.algorithms.integration.report import (
    IntegrationReport,
    ProfileModelReport,
    ProfileValidationReport,
)
from dials.algorithms.integration.validation import ValidatedMultiExpProfileModeller
from dials.algorithms.profile_model.modeller import MultiExpProfileModeller
from dials.algorithms.shoebox import MaskCode
from dials.array_family import flex
from dials.util import Sorry, phil, pprint, tabulate
from dials.util.command_line import heading
from dials.util.report import Report
from dials_algorithms_integration_integrator_ext import (
    Executor,
    JobList,
    ReflectionManager,
    max_memory_needed,
)

logger = logging.getLogger(__name__)

__all__ = [
    "Executor",
    "Integrator",
    "Integrator2D",
    "Integrator3D",
    "Integrator3DThreaded",
    "IntegratorExecutor",
    "IntegratorFlat3D",
    "IntegratorSingle2D",
    "IntegratorStills",
    "JobList",
    "Parameters",
    "Processor2D",
    "Processor3D",
    "ProcessorFlat3D",
    "ProcessorSingle2D",
    "ProcessorStills",
    "ProfileModellerExecutor",
    "ProfileValidatorExecutor",
    "ReflectionManager",
    "frame_hist",
    "generate_phil_scope",
    "hist",
    "job",
    "nframes_hist",
    "phil_scope",
]


def generate_phil_scope():
    """
    Generate the integration phil scope.

    :return: The phil scope
    """
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

        threshold = 0.95
          .type = float(value_min=0.0, value_max=1.0)
          .help = "For block size auto the block size is calculated by sorting"
                  "reflections by the number of frames they cover and then"
                  "selecting the block size to be 2*nframes[threshold] such"
                  "that 100*threshold % of reflections are guaranteed to be"
                  "fully contained in 1 block"

        force = False
          .type = bool
          .help = "If the number of processors is 1 and force is False, then the"
                  "number of blocks may be set to 1. If force is True then the"
                  "block size is always calculated."

        max_memory_usage = 0.90
          .type = float(value_min=0.0,value_max=1.0)
          .help = "The maximum percentage of available memory to use for"
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

      integrator = *auto 3d flat3d 2d single2d stills 3d_threaded
        .type = choice
        .help = "The integrator to use."
        .expert_level=3

      profile {

        fitting = True
          .type = bool
          .help = "Use profile fitting if available"

        valid_foreground_threshold = 0.75
          .type = float(value_min=0, value_max=1)
          .help = "The minimum fraction of foreground pixels that must be valid"
                  "in order for a reflection to be integrated by profile fitting."
          .expert_level = 2

        sigma_b_multiplier = 2.0
          .type = float(value_min=1.0)
          .help = "Background box expansion factor"
          .expert_level = 3

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
                  "positive value is used as the minimum permissible value."
          .type = float(value_min=0.0, value_max=1.0)

        ice_rings = False
          .help = "Set the ice ring flags"
          .type = bool
      }

      include scope dials.algorithms.integration.overlaps_filter.phil_scope

      mp {
        method = *multiprocessing drmaa sge lsf pbs
          .type = choice
          .help = "The multiprocessing method to use"

        njobs = 1
          .type = int(value_min=1)
          .help = "The number of cluster jobs to use"

        nproc = 1
          .type = int(value_min=1)
          .help = "The number of processes to use per cluster job"

        multiprocessing.n_subset_split = None
            .type = int(value_min=1)
            .help = "Number of subsets to split the reflection table for integration."
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
        [(z + 1) for b in bbox for z in range(b[4], b[5])],
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


class Parameters:
    """
    A stack of classes to represent the integration parameters
    """

    class Filter:
        """
        Filter parameters
        """

        def __init__(self):
            self.min_zeta = 0.05
            self.powder_filter = None

    class Profile:
        """
        Profile parameters
        """

        class Validation:
            def __init__(self):
                self.number_of_partitions = 2
                self.min_partition_size = 100

        def __init__(self):
            self.fitting = True
            self.sigma_b_multiplier = 2.0
            self.valid_foreground_threshold = 0.75
            self.validation = Parameters.Profile.Validation()

    def __init__(self):
        """
        Initialize
        """
        self.modelling = processor.Parameters()
        self.integration = processor.Parameters()
        self.filter = Parameters.Filter()
        self.profile = Parameters.Profile()
        self.debug_reference_filename = "reference_profiles.pickle"
        self.debug_reference_output = False

    @staticmethod
    def from_phil(params):
        """
        Convert the phil parameters
        """
        # Init the parameters
        result = Parameters()

        # Create the multi processing parameters
        mp = processor.MultiProcessing()
        mp.method = params.mp.method
        mp.nproc = params.mp.nproc
        mp.njobs = params.mp.njobs
        mp.n_subset_split = params.mp.multiprocessing.n_subset_split

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
        result.integration.integrator = params.integrator

        result.debug_reference_filename = params.debug.reference.filename
        result.debug_reference_output = params.debug.reference.output

        # Profile parameters
        result.profile.sigma_b_multiplier = params.profile.sigma_b_multiplier
        result.profile.valid_foreground_threshold = (
            params.profile.valid_foreground_threshold
        )

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


def _initialize_rotation(experiments, params, reflections):
    """
    A pre-processing class for oscillation data.
    """

    # Compute some reflection properties
    reflections.compute_zeta_multi(experiments)
    reflections.compute_d(experiments)
    reflections.compute_bbox(
        experiments, sigma_b_multiplier=params.profile.sigma_b_multiplier
    )

    # Filter the reflections by zeta
    mask = flex.abs(reflections["zeta"]) < params.filter.min_zeta
    reflections.set_flags(mask, reflections.flags.dont_integrate)

    # Filter the reflections by powder ring
    if params.filter.powder_filter is not None:
        mask = params.filter.powder_filter(reflections["d"])
        reflections.set_flags(mask, reflections.flags.in_powder_ring)


def _initialize_stills(experiments, params, reflections):
    """
    A pre-processing class for stills data.
    """

    # Compute some reflection properties
    reflections.compute_d(experiments)
    reflections.compute_bbox(
        experiments, sigma_b_multiplier=params.profile.sigma_b_multiplier
    )

    # Check the bounding boxes are all 1 frame in width
    z0, z1 = reflections["bbox"].parts()[4:6]
    assert (z1 - z0).all_eq(1), "bbox is invalid"

    # Filter the reflections by powder ring
    if params.filter.powder_filter is not None:
        mask = params.filter.powder_filter(reflections["d"])
        reflections.set_flags(mask, reflections.flags.in_powder_ring)


def _finalize(reflections, experiments, params):
    """
    A generic post-processing function.
    """
    overlaps_scope = params.integration.overlaps_filter
    if True in [
        overlaps_scope.foreground_foreground.enable,
        overlaps_scope.foreground_background.enable,
    ]:
        from dials.algorithms.integration.overlaps_filter import OverlapsFilterMultiExpt

        overlaps_filter = OverlapsFilterMultiExpt(reflections, experiments)
        if overlaps_scope.foreground_foreground.enable:
            overlaps_filter.remove_foreground_foreground_overlaps()
        if overlaps_scope.foreground_background.enable:
            overlaps_filter.remove_foreground_background_overlaps()
        reflections = overlaps_filter.refl
    return reflections, experiments


def _finalize_rotation(reflections, experiments, params):
    """
    A post-processing function for oscillation data.
    """

    reflections, experiments = _finalize(reflections, experiments, params)

    # Compute the corrections
    reflections.compute_corrections(experiments)
    return reflections, experiments


def _finalize_stills(reflections, experiments, params):
    """
    A post-processing function for stills data.
    """
    reflections, experiments = _finalize(reflections, experiments, params)

    integrated = reflections

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
        "Filtering %d reflections with at least one bad foreground pixel out of %d",
        len_all - len(integrated),
        len_all,
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
                    "Filtering %d reflections with no integrated signal (sum and variance = 0) out of %d",
                    n,
                    len(integrated),
                )
                integrated = integrated.select(integrated["intensity.sum.variance"] > 0)
            else:
                raise Sorry(
                    "Found reflections with variances == 0 but summed signal != 0"
                )

        # apply detector gain to summation variances
        integrated[
            "intensity.sum.variance"
        ] *= params.integration.summation.detector_gain
    if "background.sum.value" in integrated:
        if (integrated["background.sum.variance"] < 0).count(True) > 0:
            raise Sorry(
                "Found negative variances (background sum). Are bad pixels properly masked out?"
            )
        if (integrated["background.sum.variance"] == 0).count(True) > 0:
            logger.info(
                "Filtering %d reflections with zero background variance",
                (integrated["background.sum.variance"] == 0).count(True),
            )
            integrated = integrated.select(integrated["background.sum.variance"] > 0)
        # apply detector gain to background summation variances
        integrated[
            "background.sum.variance"
        ] *= params.integration.summation.detector_gain

    reflections = integrated

    return reflections, experiments


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
        super().__init__()

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
        logger.debug("")
        logger.debug(" Beginning modelling job %d", job.index)
        logger.info("")
        logger.info(" Frames: %d -> %d", frame0 + 1, frame1)
        logger.info("")
        logger.info(" Number of reflections")
        logger.info("  Partial:     %d", npart)
        logger.info("  Full:        %d", nfull)
        logger.info("  In ice ring: %d", nice)
        logger.info("  Total:       %d", ntot)
        logger.info("")

        # Print a histogram of reflections on frames
        if frame1 - frame0 > 1:
            logger.debug(
                " The following histogram shows the number of reflections predicted"
            )
            logger.debug(" to have all or part of their intensity on each frame.")
            logger.debug("")
            logger.debug(frame_hist(reflections["bbox"], prefix=" ", symbol="*"))
            logger.debug("")

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
        logger.debug(fmt, nmod, ntot, frame + 1)

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
        super().__init__()

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
        logger.debug("")
        logger.debug(" Beginning modelling job %d", job.index)
        logger.info("")
        logger.info(" Frames: %d -> %d", frame0, frame1)
        logger.info("")
        logger.info(" Number of reflections")
        logger.info("  Partial:     %d", npart)
        logger.info("  Full:        %d", nfull)
        logger.info("  In ice ring: %d", nice)
        logger.info("  Total:       %d", ntot)
        logger.info("")

        # Print a histogram of reflections on frames
        if frame1 - frame0 > 1:
            logger.debug(
                " The following histogram shows the number of reflections predicted"
            )
            logger.debug(" to have all or part of their intensity on each frame.")
            logger.debug("")
            logger.debug(frame_hist(reflections["bbox"], prefix=" ", symbol="*"))
            logger.debug("")

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
        logger.debug(fmt, nmod, ntot, frame + 1)

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

    def __init__(
        self, experiments, profile_fitter=None, valid_foreground_threshold=0.75
    ):
        """
        Initialize the executor

        :param experiments: The experiment list
        """
        self.experiments = experiments
        self.overlaps = None
        self.profile_fitter = profile_fitter
        self.valid_foreground_threshold = valid_foreground_threshold
        super().__init__()

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
        logger.debug("")
        logger.debug(" Beginning integration job %d", job.index)
        logger.info("")
        logger.info(" Frames: %d -> %d", frame0, frame1)
        logger.info("")
        logger.info(" Number of reflections")
        logger.info("  Partial:     %d", npart)
        logger.info("  Full:        %d", nfull)
        logger.info("  In ice ring: %d", nice)
        logger.info("  Integrate:   %d", nint)
        logger.info("  Total:       %d", ntot)
        logger.info("")

        # Print a histogram of reflections on frames
        if frame1 - frame0 > 1:
            logger.debug(
                " The following histogram shows the number of reflections predicted"
            )
            logger.debug(" to have all or part of their intensity on each frame.")
            logger.debug("")
            logger.debug(frame_hist(reflections["bbox"], prefix=" ", symbol="*"))
            logger.debug("")

        # Find any overlaps
        self.overlaps = reflections.find_overlaps(self.experiments)

    def process(self, frame, reflections):
        """
        Process the reflections on a frame

        :param frame: The frame to process
        :param reflections: The reflections to process
        """
        # Check if pixels are overloaded
        reflections.is_overloaded(self.experiments)

        # Compute the shoebox mask
        reflections.compute_mask(self.experiments)

        # Check for invalid pixels in foreground/background
        reflections.contains_invalid_pixels()

        # Exclude reflections where a high fraction of the foreground is masked
        # e.g. due to a panel edge, as this will make the fitting unreliable.
        sbox = reflections["shoebox"]
        nvalfg = sbox.count_mask_values(MaskCode.Valid | MaskCode.Foreground)
        nforeg = sbox.count_mask_values(MaskCode.Foreground)
        fraction_valid = nvalfg.as_double() / nforeg.as_double()
        selection = fraction_valid < self.valid_foreground_threshold
        reflections.set_flags(selection, reflections.flags.dont_integrate)
        logger.debug(
            f"{selection.count(True)} reflections have"
            " a fraction of valid pixels below the valid foreground threshold"
        )

        # Process the data
        reflections.compute_background(self.experiments)
        reflections.compute_centroid(self.experiments)

        reflections.compute_summed_intensity()
        if self.profile_fitter:
            reflections.compute_fitted_intensity(self.profile_fitter)

        # Compute the number of background/foreground pixels
        reflections["num_pixels.valid"] = sbox.count_mask_values(MaskCode.Valid)
        reflections["num_pixels.background"] = sbox.count_mask_values(
            MaskCode.Valid | MaskCode.Background
        )
        reflections["num_pixels.background_used"] = sbox.count_mask_values(
            MaskCode.Valid | MaskCode.Background | MaskCode.BackgroundUsed
        )
        reflections["num_pixels.foreground"] = nvalfg

        # Print some info
        fmt = " Integrated % 5d (sum) + % 5d (prf) / %5d reflections on image %d"
        nsum = reflections.get_flags(reflections.flags.integrated_sum).count(True)
        nprf = reflections.get_flags(reflections.flags.integrated_prf).count(True)
        ntot = len(reflections)
        logger.debug(fmt, nsum, nprf, ntot, frame + 1)

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


class Integrator:
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

    def fit_profiles(self):
        """Do profile fitting if appropriate.

        Sets self.profile_validation_report and self.profile_model_report.

        Returns profile_fitter (may be none)
        """
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
                profile_modellers = []
                for i in range(num_folds):
                    profile_fitter_single = MultiExpProfileModeller()  # (num_folds)
                    for expr in self.experiments:
                        profile_fitter_single.add(expr.profile.fitting_class()(expr))
                    profile_modellers.append(profile_fitter_single)

                # Create the data processor
                executor = ProfileModellerExecutor(
                    self.experiments,
                    ValidatedMultiExpProfileModeller(profile_modellers),
                )
                processor = build_processor(
                    self.ProcessorClass,
                    self.experiments,
                    reference,
                    self.params.modelling,
                )
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

                # Dump reference profiles
                if self.params.debug_reference_output:
                    reference_debug = []
                    for i in range(len(finalized_profile_fitter)):
                        m = finalized_profile_fitter[i]
                        p = []
                        for j in range(len(m)):
                            try:
                                p.append(
                                    {
                                        "data": m.data(j),
                                        "mask": m.mask(j),
                                        "coord": m.coord(j),
                                        "n_reflections": m.n_reflections(j),
                                    }
                                )
                            except Exception:
                                p.append(None)
                        reference_debug.append(p)
                    with open(self.params.debug_reference_filename, "wb") as outfile:
                        pickle.dump(reference_debug, outfile)

                # Print profiles
                for i in range(len(finalized_profile_fitter)):
                    m = finalized_profile_fitter[i]
                    logger.debug("")
                    logger.debug("Profiles for experiment %d", i)
                    for j in range(len(m)):
                        logger.debug("Profile %d", j)
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
                    processor = build_processor(
                        self.ProcessorClass,
                        self.experiments,
                        reference,
                        self.params.modelling,
                    )
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
        return profile_fitter

    def integrate(self):
        """
        Integrate the data
        """
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

        # Print the summary
        logger.info(
            " Processing the following experiments:\n"
            "\n"
            " Experiments: %d\n"
            " Beams:       %d\n"
            " Detectors:   %d\n"
            " Goniometers: %d\n"
            " Scans:       %d\n"
            " Crystals:    %d\n"
            " Imagesets:   %d\n",
            len(self.experiments),
            len(self.experiments.beams()),
            len(self.experiments.detectors()),
            len(self.experiments.goniometers()),
            len(self.experiments.scans()),
            len(self.experiments.crystals()),
            len(self.experiments.imagesets()),
        )

        # Initialize the reflections
        self.initialize_reflections(self.experiments, self.params, self.reflections)

        # Check if we want to do some profile fitting
        profile_fitter = self.fit_profiles()

        logger.info("=" * 80)
        logger.info("")
        logger.info(heading("Integrating reflections"))
        logger.info("")

        # Create the data processor
        executor = IntegratorExecutor(
            self.experiments,
            profile_fitter,
            self.params.profile.valid_foreground_threshold,
        )

        # determine the max memory needed during integration
        def _determine_max_memory_needed(experiments, reflections):
            max_needed = 0
            for imageset in experiments.imagesets():
                # find all experiments belonging to that imageset, as each
                # imageset is processed as a whole for integration.
                if all(experiments.identifiers()):
                    expt_ids = [
                        experiment.identifier
                        for experiment in experiments
                        if experiment.imageset == imageset
                    ]
                    subset = reflections.select_on_experiment_identifiers(expt_ids)
                else:
                    subset = flex.reflection_table()
                    for j, experiment in enumerate(experiments):
                        if experiment.imageset == imageset:
                            subset.extend(reflections.select(reflections["id"] == j))
                try:
                    if imageset.get_scan():
                        frame0, frame1 = imageset.get_scan().get_array_range()
                    else:
                        raise RuntimeError
                except RuntimeError:  # catch DXTBX_ASSERT if no scan in imageset
                    frame0, frame1 = (0, len(imageset))
                flatten = self.params.integration.integrator == "flat3d"
                max_needed = max(
                    max_memory_needed(subset, frame0, frame1, flatten),
                    max_needed,
                )
            assert max_needed > 0, "Could not determine memory requirements"
            return max_needed

        def _iterative_table_split(tables, experiments, available_memory):
            split_tables = []
            for table in tables:
                mem_needed = _determine_max_memory_needed(experiments, table)
                if mem_needed > available_memory:
                    n_to_split = int(math.ceil(mem_needed / available_memory))
                    flex.set_random_seed(0)
                    split_tables.extend(table.random_split(n_to_split))
                else:
                    split_tables.append(table)
            if len(split_tables) == len(tables):
                # nothing was split, all passed memory check
                return split_tables
            # some tables were split - so need to check again that all are ok
            return _iterative_table_split(split_tables, experiments, available_memory)

        def _run_processor(reflections):
            processor = build_processor(
                self.ProcessorClass,
                self.experiments,
                reflections,
                self.params.integration,
            )
            processor.executor = executor
            # Process the reflections
            reflections, _, time_info = processor.process()
            return reflections, time_info

        if self.params.integration.mp.method != "multiprocessing":
            self.reflections, time_info = _run_processor(self.reflections)
        else:
            # need to do a memory check and decide whether to split table
            available_immediate, _, __ = assess_available_memory(
                self.params.integration
            )

            #  here don't consider nproc as the processor will reduce nproc to 1
            # if necessary, only want to split if we can't even process with
            # nproc = 1

            if self.params.integration.mp.n_subset_split:
                tables = self.reflections.random_split(
                    self.params.integration.mp.n_subset_split
                )
            else:
                tables = _iterative_table_split(
                    [self.reflections],
                    self.experiments,
                    available_immediate,
                )

            if len(tables) == 1:
                # will not fail a memory check in the processor, so proceed
                self.reflections, time_info = _run_processor(self.reflections)
            else:
                # Split the reflections and process by performing multiple
                # passes over each imageset
                time_info = TimingInfo()
                reflections = flex.reflection_table()

                logger.info(
                    """Predicted maximum memory needed exceeds available memory.
Splitting reflection table into %s subsets for processing
""",
                    len(tables),
                )
                for i, table in enumerate(tables):
                    logger.info("Processing subset %s of reflection table", i + 1)
                    processed, this_time_info = _run_processor(table)
                    reflections.extend(processed)
                    time_info += this_time_info
                self.reflections = reflections

        # Finalize the reflections
        self.reflections, self.experiments = self.finalize_reflections(
            self.reflections, self.experiments, self.params
        )

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
        result = Report()
        if self.profile_model_report is not None:
            result.combine(self.profile_model_report)
        result.combine(self.integration_report)
        return result

    def summary(self, block_size, block_size_units):
        """Print a summary of the integration stuff."""
        # Compute the task table
        if self._experiments.all_stills():
            rows = [["#", "Group", "Frame From", "Frame To"]]
            for i in range(len(self)):
                job = self._manager.job(i)
                group = job.index()
                f0, f1 = job.frames()
                rows.append([str(i), str(group), str(f0), str(f1)])
        elif self._experiments.all_sequences():
            rows = [["#", "Group", "Frame From", "Frame To", "Angle From", "Angle To"]]
            for i in range(len(self)):
                job = self._manager.job(i)
                group = job.index()
                expr = job.expr()
                f0, f1 = job.frames()
                scan = self._experiments[expr[0]].scan
                p0 = scan.get_angle_from_array_index(f0)
                p1 = scan.get_angle_from_array_index(f1)
                rows.append(
                    [str(i), str(group), str(f0 + 1), str(f1), str(p0), str(p1)]
                )
        else:
            raise RuntimeError("Experiments must be all sequences or all stills")
        return tabulate(rows, headers="firstrow")


class Integrator3D(Integrator):
    """
    Integrator for 3D algorithms
    """

    initialize_reflections = staticmethod(_initialize_rotation)
    ProcessorClass = Processor3D
    finalize_reflections = staticmethod(_finalize_rotation)


class IntegratorFlat3D(Integrator):
    """
    Integrator for flattened 3D algorithms
    """

    initialize_reflections = staticmethod(_initialize_rotation)
    ProcessorClass = ProcessorFlat3D
    finalize_reflections = staticmethod(_finalize_rotation)


class Integrator2D(Integrator):
    """
    Integrator for 2D algorithms
    """

    initialize_reflections = staticmethod(_initialize_rotation)
    ProcessorClass = Processor2D
    finalize_reflections = staticmethod(_finalize_rotation)


class IntegratorSingle2D(Integrator):
    """
    Integrator for 2D algorithms on a single image
    """

    initialize_reflections = staticmethod(_initialize_rotation)
    ProcessorClass = ProcessorSingle2D
    finalize_reflections = staticmethod(_finalize_rotation)


class IntegratorStills(Integrator):
    """
    Integrator for still algorithms
    """

    initialize_reflections = staticmethod(_initialize_stills)
    ProcessorClass = ProcessorStills
    finalize_reflections = staticmethod(_finalize_stills)


class Integrator3DThreaded:
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
        self.reflections.compute_bbox(
            self.experiments,
            sigma_b_multiplier=self.params.integration.profile.sigma_b_multiplier,
        )

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
            fmt,
            len(self.experiments),
            len(self.experiments.beams()),
            len(self.experiments.detectors()),
            len(self.experiments.goniometers()),
            len(self.experiments.scans()),
            len(self.experiments.crystals()),
            len(self.experiments.imagesets()),
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
        result = Report()
        if self.profile_model_report is not None:
            result.combine(self.profile_model_report)
        result.combine(self.integration_report)
        return result

    def summary(self, block_size, block_size_units):
        """Print a summary of the integration stuff."""
        # Compute the task table
        if self._experiments.all_stills():
            rows = [["#", "Group", "Frame From", "Frame To"]]
            for i in range(len(self)):
                job = self._manager.job(i)
                group = job.index()
                f0, f1 = job.frames()
                rows.append([str(i), str(group), str(f0), str(f1)])
        elif self._experiments.all_sequences():
            rows = [["#", "Group", "Frame From", "Frame To", "Angle From", "Angle To"]]
            for i in range(len(self)):
                job = self._manager.job(i)
                group = job.index()
                expr = job.expr()
                f0, f1 = job.frames()
                scan = self._experiments[expr[0]].scan
                p0 = scan.get_angle_from_array_index(f0)
                p1 = scan.get_angle_from_array_index(f1)
                rows.append(
                    [str(i), str(group), str(f0 + 1), str(f1), str(p0), str(p1)]
                )
        else:
            raise RuntimeError("Experiments must be all sequences or all stills")
        return tabulate(rows, headers="firstrow")


def create_integrator(params, experiments, reflections):
    """
    Create an integrator object with a given configuration.

    :param params: The input phil parameters
    :param experiments: The list of experiments
    :param reflections: The reflections to integrate
    :return: An integrator object
    """
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
    if params.integration.lookup.mask and isinstance(
        params.integration.lookup.mask, str
    ):
        with open(params.integration.lookup.mask, "rb") as infile:
            params.integration.lookup.mask = pickle.load(infile, encoding="bytes")

    # Set algorithms as reflection table defaults
    BackgroundAlgorithm = dials.extensions.Background.load(
        params.integration.background.algorithm
    )
    flex.reflection_table.background_algorithm = functools.partial(
        BackgroundAlgorithm, params
    )
    CentroidAlgorithm = dials.extensions.Centroid.load(
        params.integration.centroid.algorithm
    )
    flex.reflection_table.centroid_algorithm = functools.partial(
        CentroidAlgorithm, params
    )

    # Get the classes we need
    if params.integration.integrator == "auto":
        if experiments.all_stills():
            params.integration.integrator = "stills"
        else:
            params.integration.integrator = "3d"
    IntegratorClass = {
        "3d": Integrator3D,
        "flat3d": IntegratorFlat3D,
        "2d": Integrator2D,
        "single2d": IntegratorSingle2D,
        "stills": IntegratorStills,
        "3d_threaded": Integrator3DThreaded,
    }.get(params.integration.integrator)
    if not IntegratorClass:
        raise ValueError(f"Unknown integration type {params.integration.integrator}")

    # Remove scan if stills
    if experiments.all_stills():
        for experiment in experiments:
            experiment.scan = None

    # Return an instantiation of the class
    return IntegratorClass(experiments, reflections, params)
