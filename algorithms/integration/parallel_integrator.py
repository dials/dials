from __future__ import annotations

import logging
import math
import pickle

import psutil

from libtbx import Auto

import dials.algorithms.integration
from dials.algorithms.integration.processor import NullTask, execute_parallel_task
from dials.array_family import flex
from dials.util import tabulate
from dials.util.mp import multi_node_parallel_map

# Need this import first because loads extension that parallel_integrator_ext
# relies on - it assumes the binding for EmpiricalProfileModeller exists
import dials.algorithms.profile_model.modeller  # noqa: F401 # isort: split

from dials_algorithms_integration_parallel_integrator_ext import (
    GaussianRSIntensityCalculator,
    GaussianRSMaskCalculator,
    GaussianRSMultiCrystalMaskCalculator,
    GaussianRSMultiCrystalReferenceProfileData,
    GaussianRSReferenceCalculator,
    GaussianRSReferenceProfileData,
    GLMBackgroundCalculator,
    Logger,
    MultiThreadedIntegrator,
    MultiThreadedReferenceProfiler,
    ReferenceProfileData,
    SimpleBackgroundCalculator,
    SimpleBlockList,
    SimpleReflectionManager,
)

__all__ = [
    "BackgroundCalculatorFactory",
    "GLMBackgroundCalculator",
    "GaussianRSIntensityCalculator",
    "GaussianRSMaskCalculator",
    "GaussianRSMultiCrystalMaskCalculator",
    "GaussianRSMultiCrystalReferenceProfileData",
    "GaussianRSReferenceCalculator",
    "GaussianRSReferenceProfileData",
    "IntegrationJob",
    "IntegrationManager",
    "IntegratorProcessor",
    "IntensityCalculatorFactory",
    "Logger",
    "MaskCalculatorFactory",
    "MultiThreadedIntegrator",
    "MultiThreadedReferenceProfiler",
    "ReferenceCalculatorFactory",
    "ReferenceCalculatorJob",
    "ReferenceCalculatorManager",
    "ReferenceCalculatorProcessor",
    "ReferenceProfileData",
    "SimpleBackgroundCalculator",
    "SimpleBlockList",
    "SimpleReflectionManager",
]

logger = logging.getLogger(__name__)


class MaskCalculatorFactory:
    """
    A factory function to return a mask calculator object
    """

    @staticmethod
    def create(experiments, params=None):
        """
        Select the mask calculator
        """
        from dials.algorithms.profile_model.gaussian_rs.algorithm import (
            GaussianRSMaskCalculatorFactory,
        )

        # Get the parameters
        if params is None:
            from dials.command_line.integrate import phil_scope

            params = phil_scope.extract()

        # Select the factory function
        selection = params.profile.algorithm
        if selection == "gaussian_rs":
            algorithm = GaussianRSMaskCalculatorFactory.create(experiments)
        else:
            raise RuntimeError("Unknown profile model algorithm")

        # Create the mask algorithm
        return algorithm


class BackgroundCalculatorFactory:
    """
    A factory function to return a background calculator object
    """

    @staticmethod
    def create(experiments, params=None):
        """
        Select the background calculator
        """
        from dials.algorithms.background.glm.algorithm import (
            GLMBackgroundCalculatorFactory,
        )
        from dials.algorithms.background.gmodel.algorithm import (
            GModelBackgroundCalculatorFactory,
        )
        from dials.algorithms.background.simple.algorithm import (
            SimpleBackgroundCalculatorFactory,
        )

        # Get the parameters
        if params is None:
            from dials.command_line.integrate import phil_scope

            params = phil_scope.extract()

        # Select the factory function
        selection = params.integration.background.algorithm
        if selection == "simple":

            # Get parameters
            params = params.integration.background.simple

            # Create some keyword parameters
            kwargs = {
                "model": params.model.algorithm,
                "outlier": params.outlier.algorithm,
                "min_pixels": params.min_pixels,
            }

            # Create all the keyword parameters
            if params.outlier.algorithm == "null":
                pass
            elif params.outlier.algorithm == "truncated":
                kwargs["lower"] = params.outlier.truncated.lower
                kwargs["upper"] = params.outlier.truncated.upper
            elif params.outlier.algorithm == "nsigma":
                kwargs["lower"] = params.outlier.nsigma.lower
                kwargs["upper"] = params.outlier.nsigma.upper
            elif params.outlier.algorithm == "normal":
                kwargs["min_pixels"] = params.outlier.normal.min_pixels
            elif params.outlier.algorithm == "plane":
                kwargs["fraction"] = params.outlier.plane.fraction
                kwargs["n_sigma"] = params.outlier.plane.n_sigma
            elif params.outlier.algorithm == "tukey":
                kwargs["lower"] = params.outlier.tukey.lower
                kwargs["upper"] = params.outlier.tukey.upper

            # Create the algorithm
            algorithm = SimpleBackgroundCalculatorFactory.create(experiments, **kwargs)

        elif selection == "glm":

            # Get the parameters
            params = params.integration.background.glm

            # Create the algorithm
            algorithm = GLMBackgroundCalculatorFactory.create(
                experiments,
                model=params.model.algorithm,
                tuning_constant=params.robust.tuning_constant,
                min_pixels=params.min_pixels,
            )

        elif selection == "gmodel":

            # Get the parameters
            params = params.integration.background.gmodel

            # Create the algorithm
            algorithm = GModelBackgroundCalculatorFactory.create(
                experiments,
                model=params.model,
                robust=params.robust.algorithm,
                min_pixels=params.min_pixels,
            )

        else:
            raise RuntimeError("Unknown background algorithm")

        # Return the background calculator
        return algorithm


class IntensityCalculatorFactory:
    """
    A factory function to return an intensity calculator object
    """

    @staticmethod
    def create(experiments, reference_profiles, params=None):
        """
        Select the intensity calculator
        """
        from dials.algorithms.profile_model.gaussian_rs.algorithm import (
            GaussianRSIntensityCalculatorFactory,
        )

        # Get the parameters
        if params is None:
            from dials.command_line.integrate import phil_scope

            params = phil_scope.extract()

        # Select the factory function
        selection = params.profile.algorithm
        if selection == "gaussian_rs":

            # Get the parameters
            params = params.profile.gaussian_rs.fitting

            # Check for detector space
            if params.fit_method == "reciprocal_space":
                detector_space = False
            elif params.fit_method == "detector_space":
                detector_space = True
            else:
                raise RuntimeError(f"Unknown fit method: {params.fit_method}")

            # Create the algorithm
            algorithm = GaussianRSIntensityCalculatorFactory.create(
                reference_profiles,
                detector_space=detector_space,
                deconvolution=params.detector_space.deconvolution,
            )

        else:
            raise RuntimeError("Unknown profile model algorithm")

        # Return the algorithm
        return algorithm


class ReferenceCalculatorFactory:
    """
    A factory function to return an reference calculator object
    """

    @staticmethod
    def create(experiments, params=None):
        """
        Select the reference calculator
        """
        from dials.algorithms.profile_model.gaussian_rs.algorithm import (
            GaussianRSReferenceCalculatorFactory,
        )

        # Get the parameters
        if params is None:
            from dials.command_line.integrate import phil_scope

            params = phil_scope.extract()

        # Select the factory function
        selection = params.profile.algorithm
        if selection == "gaussian_rs":

            # Get the parameters
            params = params.profile.gaussian_rs.fitting

            # Create the algorithm
            algorithm = GaussianRSReferenceCalculatorFactory.create(
                experiments,
                grid_size=params.grid_size,
                scan_step=params.scan_step,
                grid_method=params.grid_method,
            )

        else:
            raise RuntimeError("Unknown profile model algorithm")

        # Return the algorithm
        return algorithm


def _assert_enough_memory(required_memory, max_memory_usage):
    """
    Check there is enough memory available or fail

    :param required_memory: The required number of bytes
    :param max_memory_usage: The maximum memory usage allowed
    """
    total_memory = psutil.virtual_memory().total
    assert max_memory_usage > 0.0, "maximum memory usage must be > 0"
    assert max_memory_usage <= 1.0, "maximum memory usage must be <= 1"
    limit_memory = total_memory * max_memory_usage
    if required_memory > limit_memory:
        raise RuntimeError(
            f"""
    There was a problem allocating memory for image data. Possible solutions
    include increasing the percentage of memory allowed for shoeboxes or
    decreasing the block size. This could also be caused by a highly mosaic
    crystal model - is your crystal really this mosaic?
      Total system memory: {total_memory / 1000000000.0:.1f} GB
      Limit image memory: {limit_memory / 1000000000.0:.1f} GB
      Required image memory: {required_memory / 1000000000.0:.1f} GB
    """
        )
    else:
        logger.info("Allocating %.1f MB memory", required_memory / 1e6)


class IntegrationJob:
    """
    A class to represent an integration job
    """

    def __init__(self, index, job, experiments, reflections, reference, params=None):
        """
        Initialise the task.

        :param index: The index of the processing job
        :param experiments: The list of experiments
        :param reflections: The list of reflections
        :param params: The processing parameters
        :param job: The frames to integrate
        :param flatten: Flatten the shoeboxes
        :param executor: The executor class
        """

        # Get the parameters
        if params is None:
            from dials.command_line.integrate import phil_scope

            params = phil_scope.extract()

        assert len(reflections) > 0, "Zero reflections given"
        assert (
            params.integration.block.max_memory_usage > 0.0
        ), "Max memory % must be > 0"
        assert (
            params.integration.block.max_memory_usage <= 1.0
        ), "Max memory % must be < 1"
        self.index = index
        self.job = job
        self.experiments = experiments
        self.reflections = reflections
        self.reference = reference
        self.params = params

    def __call__(self):
        """
        Do the processing.

        :return: The processed data
        """
        from dials.algorithms.integration.processor import job

        # Set the global process ID
        job.index = self.index

        # Check all reflections have same imageset and get it
        imageset = self.experiments[0].imageset
        if not all(e.imageset == imageset for e in self.experiments):
            raise RuntimeError("Task can only handle 1 imageset")

        # Get the sub imageset
        frame0, frame1 = self.job
        try:
            frame10, frame11 = imageset.get_array_range()
        except Exception:
            frame10, frame11 = (0, len(imageset))
        try:
            assert frame0 < frame1
            assert frame10 < frame11
            assert frame0 >= frame10
            assert frame1 <= frame11
            index0 = frame0 - frame10
            index1 = index0 + (frame1 - frame0)
            assert index0 < index1
            assert index0 >= 0
            assert index1 <= len(imageset)
            imageset = imageset[index0:index1]
        except Exception:
            raise RuntimeError("Programmer Error: bad array range")

        # Check the memory requirements
        _assert_enough_memory(
            self.compute_required_memory(imageset),
            self.params.integration.block.max_memory_usage,
        )

        # Integrate
        self.integrate(imageset)

        # Write some debug files
        self.write_debug_files()

        # Return the result
        return dials.algorithms.integration.Result(
            index=self.index,
            reflections=self.reflections,
            data=None,
            read_time=0,
            extract_time=0,
            process_time=0,
            total_time=0,
        )

    def compute_required_memory(self, imageset):
        """
        Compute the required memory
        """
        return MultiThreadedIntegrator.compute_required_memory(
            imageset, self.params.integration.block.size
        )

    def integrate(self, imageset):
        """
        Integrate the reflections
        """
        from dials.algorithms.integration.integrator import frame_hist

        # Compute the partiality
        self.reflections.compute_partiality(self.experiments)

        # Get some info
        EPS = 1e-7
        full_value = 0.997300203937 - EPS
        fully_recorded = self.reflections["partiality"] > full_value
        npart = fully_recorded.count(False)
        nfull = fully_recorded.count(True)
        select_ice = self.reflections.get_flags(self.reflections.flags.in_powder_ring)
        select_int = ~self.reflections.get_flags(self.reflections.flags.dont_integrate)
        nice = select_ice.count(True)
        nint = select_int.count(True)
        ntot = len(self.reflections)
        frame0, frame1 = imageset.get_scan().get_array_range()

        # Write some output
        logger.info(" Beginning integration job %d", self.index)
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
            logger.info(
                " The following histogram shows the number of reflections predicted"
            )
            logger.info(" to have all or part of their intensity on each frame.")
            logger.info("")
            logger.info(
                frame_hist(
                    self.reflections["bbox"].select(select_int), prefix=" ", symbol="*"
                )
            )
            logger.info("")

        # Construct the mask algorithm
        compute_mask = MaskCalculatorFactory.create(self.experiments, self.params)

        # Construct the background algorithm
        compute_background = BackgroundCalculatorFactory.create(
            self.experiments, self.params
        )

        # Construct the intensity algorithm
        compute_intensity = IntensityCalculatorFactory.create(
            self.experiments, self.reference, self.params
        )

        # Call the multi threaded integrator
        integrator = MultiThreadedIntegrator(
            reflections=self.reflections,
            imageset=imageset,
            compute_mask=compute_mask,
            compute_background=compute_background,
            compute_intensity=compute_intensity,
            logger=Logger(logger),
            nthreads=self.params.integration.mp.nproc,
            buffer_size=self.params.integration.block.size,
            use_dynamic_mask=self.params.integration.use_dynamic_mask,
            debug=self.params.integration.debug.output,
        )

        # Assign the reflections
        self.reflections = integrator.reflections()

    def write_debug_files(self):
        """
        Write some debug output
        """

        # Optionally save the shoeboxes
        debug = self.params.integration.debug
        if debug.output and debug.separate_files:
            output = self.reflections
            if debug.select is not None:
                output = output.select(debug.select(output))
            if debug.split_experiments:
                output = output.split_by_experiment_id()
                for table in output:
                    i = table["id"][0]
                    table.as_file("shoeboxes_%d_%d.refl" % (self.index, i))
            else:
                output.as_file("shoeboxes_%d.refl" % self.index)

        # Delete the shoeboxes
        if debug.separate_files or not debug.output:
            del self.reflections["shoebox"]


class IntegrationManager:
    """
    A class to manage processing book-keeping
    """

    def __init__(self, experiments, reflections, reference, params):
        """
        Initialise the manager.

        :param experiments: The list of experiments
        :param reflections: The list of reflections
        :param reference: The reference profiles
        :param params: The phil parameters
        """

        # Save some data
        self.experiments = experiments
        self.reflections = reflections
        self.reference = reference

        # Save some parameters
        self.params = params

        # Set the finalized flag to False
        self.finalized = False

        # Initialise the timing information
        # self.time = TimingInfo()

        self.initialize()

    def initialize(self):
        """
        Initialise the processing
        """
        # Ensure the reflections contain bounding boxes
        assert "bbox" in self.reflections, "Reflections have no bbox"

        # Compute the block size and jobs
        self.compute_blocks()
        self.compute_jobs()
        self.reflections = split_partials_over_boundaries(
            self.reflections, self.params.integration.block.size
        )

        # Create the reflection manager
        self.manager = SimpleReflectionManager(
            self.blocks, self.reflections, self.params.integration.mp.njobs
        )

    def task(self, index):
        """
        Get a task.
        """
        frames = self.manager.job(index)
        experiments = self.experiments
        reference = self.reference
        reflections = self.manager.split(index)
        if len(reflections) == 0:
            logger.warning("No reflections in job %d ***", index)
            task = NullTask(index=index, reflections=reflections)
        else:
            task = IntegrationJob(
                index=index,
                job=frames,
                experiments=experiments,
                reflections=reflections,
                reference=reference,
                params=self.params,
            )
        return task

    def tasks(self):
        """
        Iterate through the tasks.
        """
        for i in range(len(self)):
            yield self.task(i)

    def accumulate(self, result):
        """Accumulate the results."""
        self.manager.accumulate(result.index, result.reflections)
        # self.time.read += result.read_time
        # self.time.extract += result.extract_time
        # self.time.process += result.process_time
        # self.time.total += result.total_time

    def finalize(self):
        """
        Finalize the processing and finish.
        """
        # Check manager is finished
        assert self.manager.finished(), "Manager is not finished"

        self.finalized = True

    def result(self):
        """
        Return the result.

        :return: The result
        """
        assert self.finalized, "Manager is not finalized"
        return self.manager.data()

    def finished(self):
        """
        Return if all tasks have finished.

        :return: True/False all tasks have finished
        """
        return self.finalized and self.manager.finished()

    def __len__(self):
        """
        Return the number of tasks.

        :return: the number of tasks
        """
        return len(self.manager)

    def compute_max_block_size(self):
        """
        Compute the required memory
        """
        total_memory = psutil.virtual_memory().available
        max_memory_usage = self.params.integration.block.max_memory_usage
        assert max_memory_usage > 0.0, "maximum memory usage must be > 0"
        assert max_memory_usage <= 1.0, "maximum memory usage must be <= 1"
        limit_memory = int(math.floor(total_memory * max_memory_usage))
        return MultiThreadedIntegrator.compute_max_block_size(
            self.experiments[0].imageset, max_memory_usage=limit_memory
        )

    def compute_blocks(self):
        """
        Compute the processing block size.
        """
        block = self.params.integration.block
        max_block_size = self.compute_max_block_size()
        if block.size in [Auto, "auto", "Auto"]:
            assert block.threshold > 0, "Threshold must be > 0"
            assert block.threshold <= 1.0, "Threshold must be < 1"
            nframes = sorted([b[5] - b[4] for b in self.reflections["bbox"]])
            cutoff = int(block.threshold * len(nframes))
            block_size = nframes[cutoff]
            if block_size > max_block_size:
                logger.warning(
                    "Computed block size (%s) > maximum block size (%s).",
                    block_size,
                    max_block_size,
                )
                logger.warning(
                    "Setting block size to maximum; some reflections may be partial"
                )
                block_size = max_block_size
        else:
            scan = self.experiments[0].scan
            if block.units == "radians":
                phi0, dphi = scan.get_oscillation(deg=False)
                block_size = int(math.ceil(block.size / dphi))
            elif block.units == "degrees":
                phi0, dphi = scan.get_oscillation()
                block_size = int(math.ceil(block.size / dphi))
            elif block.units == "frames":
                block_size = int(math.ceil(block.size))
            else:
                raise RuntimeError(f"Unknown block_size unit {block.units!r}")
            if block_size > max_block_size:
                raise RuntimeError(
                    """
          The requested block size (%s) is larger than the maximum allowable block
          size (%s). Either decrease the requested block size or increase the
          amount of available memory.
        """
                    % (block_size, max_block_size)
                )
        block.size = block_size
        block.units = "frames"

    def compute_jobs(self):
        """
        Compute the jobs
        """
        imageset = self.experiments[0].imageset
        array_range = imageset.get_array_range()
        block = self.params.integration.block
        assert block.units == "frames"
        assert block.size > 0
        self.blocks = SimpleBlockList(array_range, block.size)
        assert len(self.blocks) > 0, "Invalid number of jobs"

    def summary(self):
        """
        Get a summary of the processing
        """
        # Compute the task table
        if self.experiments.all_stills():
            rows = [["#", "Group", "Frame From", "Frame To", "# Reflections"]]
            for i in range(len(self)):
                job = self.manager.job(i)
                group = job.index()
                f0, f1 = job.frames()
                n = self.manager.num_reflections(i)
                rows.append([str(i), str(group), str(f0), str(f1), str(n)])
        elif self.experiments.all_sequences():
            rows = [
                [
                    "#",
                    "Frame From",
                    "Frame To",
                    "Angle From",
                    "Angle To",
                    "# Reflections",
                ]
            ]
            for i in range(len(self)):
                f0, f1 = self.manager.job(i)
                scan = self.experiments[0].scan
                p0 = scan.get_angle_from_array_index(f0)
                p1 = scan.get_angle_from_array_index(f1)
                n = self.manager.num_reflections(i)
                rows.append([str(i), str(f0), str(f1), str(p0), str(p1), str(n)])
        else:
            raise RuntimeError("Experiments must be all sequences or all stills")

        # The job table
        task_table = tabulate(rows, headers="firstrow")

        # The format string
        if self.params.integration.block.size is None:
            block_size = "auto"
        else:
            block_size = str(self.params.integration.block.size)
        fmt = (
            "Processing reflections in the following blocks of images:\n"
            "\n"
            " block_size: %s %s\n"
            "\n"
            "%s\n"
        )
        return fmt % (block_size, self.params.integration.block.units, task_table)


class ReferenceCalculatorJob:
    """
    A class to represent an integration job
    """

    def __init__(self, index, job, experiments, reflections, params=None):
        """
        Initialise the task.

        :param index: The index of the processing job
        :param experiments: The list of experiments
        :param reflections: The list of reflections
        :param params: The processing parameters
        :param job: The frames to integrate
        :param flatten: Flatten the shoeboxes
        :param executor: The executor class
        """

        # Get the parameters
        if params is None:
            from dials.command_line.integrate import phil_scope

            params = phil_scope.extract()

        assert len(reflections) > 0, "Zero reflections given"
        assert (
            params.integration.block.max_memory_usage > 0.0
        ), "Max memory % must be > 0"
        assert (
            params.integration.block.max_memory_usage <= 1.0
        ), "Max memory % must be < 1"
        self.index = index
        self.job = job
        self.experiments = experiments
        self.reflections = reflections
        self.params = params

    def __call__(self):
        """
        Do the processing.

        :return: The processed data
        """
        from dials.algorithms.integration.processor import job

        # Set the global process ID
        job.index = self.index

        # Check all reflections have same imageset and get it
        imageset = self.experiments[0].imageset
        if not all(e.imageset == imageset for e in self.experiments):
            raise RuntimeError("Task can only handle 1 imageset")

        # Get the sub imageset
        frame0, frame1 = self.job
        try:
            frame10, frame11 = imageset.get_array_range()
        except Exception:
            frame10, frame11 = (0, len(imageset))
        try:
            assert frame0 < frame1
            assert frame10 < frame11
            assert frame0 >= frame10
            assert frame1 <= frame11
            index0 = frame0 - frame10
            index1 = index0 + (frame1 - frame0)
            assert index0 < index1
            assert index0 >= 0
            assert index1 <= len(imageset)
            imageset = imageset[index0:index1]
        except Exception:
            raise RuntimeError("Programmer Error: bad array range")

        # Check the memory requirements
        _assert_enough_memory(
            self.compute_required_memory(imageset),
            self.params.integration.block.max_memory_usage,
        )

        # Integrate
        self.compute_reference_profiles(imageset)

        # Write some debug files
        self.write_debug_files()

        # Return the result
        return dials.algorithms.integration.Result(
            index=self.index,
            reflections=self.reflections,
            data=self.reference,
            read_time=0,
            extract_time=0,
            process_time=0,
            total_time=0,
        )

    def compute_required_memory(self, imageset):
        return MultiThreadedIntegrator.compute_required_memory(
            imageset, self.params.integration.block.size
        )

    def compute_reference_profiles(self, imageset):
        """
        Integrate the reflections
        """
        from dials.algorithms.integration.integrator import frame_hist

        # Compute the partiality
        self.reflections.compute_partiality(self.experiments)

        # Get some info
        EPS = 1e-7
        full_value = 0.997300203937 - EPS
        fully_recorded = self.reflections["partiality"] > full_value
        npart = fully_recorded.count(False)
        nfull = fully_recorded.count(True)
        select_ice = self.reflections.get_flags(self.reflections.flags.in_powder_ring)
        select_int = ~self.reflections.get_flags(self.reflections.flags.dont_integrate)
        nice = select_ice.count(True)
        nint = select_int.count(True)
        ntot = len(self.reflections)
        frame0, frame1 = imageset.get_scan().get_array_range()

        # Write some output
        logger.info(" Beginning integration job %d", self.index)
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
            logger.info(
                " The following histogram shows the number of reflections predicted"
            )
            logger.info(" to have all or part of their intensity on each frame.")
            logger.info("")
            logger.info(
                frame_hist(
                    self.reflections["bbox"].select(select_int), prefix=" ", symbol="*"
                )
            )
            logger.info("")

        # Construct the mask algorithm
        compute_mask = MaskCalculatorFactory.create(self.experiments, self.params)

        # Construct the background algorithm
        compute_background = BackgroundCalculatorFactory.create(
            self.experiments, self.params
        )

        # Construct the intensity algorithm
        compute_reference = ReferenceCalculatorFactory.create(
            self.experiments, self.params
        )

        # Call the multi threaded integrator
        reference_calculator = MultiThreadedReferenceProfiler(
            reflections=self.reflections,
            imageset=imageset,
            compute_mask=compute_mask,
            compute_background=compute_background,
            compute_reference=compute_reference,
            logger=Logger(logger),
            nthreads=self.params.integration.mp.nproc,
            buffer_size=self.params.integration.block.size,
            use_dynamic_mask=self.params.integration.use_dynamic_mask,
            debug=self.params.integration.debug.output,
        )

        # Assign the reflections
        self.reflections = reference_calculator.reflections()

        # Assign the reference profiles
        self.reference = compute_reference

        # Write some log output
        fmt = "Used %d / %d reflections to create reference profiles"
        dont_integrate = self.reflections.get_flags(
            self.reflections.flags.dont_integrate
        )
        used_in_modelling = self.reflections.get_flags(
            self.reflections.flags.used_in_modelling
        )
        n_tot = dont_integrate.count(False)
        n_mod = (used_in_modelling & ~dont_integrate).count(True)
        logger.info("")
        logger.info(fmt, n_mod, n_tot)

    def write_debug_files(self):
        """
        Write some debug output
        """

        # Optionally save the shoeboxes
        debug = self.params.integration.debug
        if debug.output and debug.separate_files:
            output = self.reflections
            if debug.select is not None:
                output = output.select(debug.select(output))
            if debug.split_experiments:
                output = output.split_by_experiment_id()
                for table in output:
                    i = table["id"][0]
                    table.as_file("shoeboxes_%d_%d.refl" % (self.index, i))
            else:
                output.as_file("shoeboxes_%d.refl" % self.index)

        # Delete the shoeboxes
        if debug.separate_files or not debug.output:
            del self.reflections["shoebox"]


class ReferenceCalculatorManager:
    """
    A class to manage processing book-keeping
    """

    def __init__(self, experiments, reflections, params):
        """
        Initialise the manager.

        :param experiments: The list of experiments
        :param reflections: The list of reflections
        :param params: The phil parameters
        """

        # Save some data
        self.experiments = experiments
        self.reflections = reflections
        self.reference = None

        # Save some parameters
        self.params = params

        # Set the finalized flag to False
        self.finalized = False

        # Initialise the timing information
        # self.time = TimingInfo()

        self.initialize()

    def initialize(self):
        """
        Initialise the processing
        """
        # Ensure the reflections contain bounding boxes
        assert "bbox" in self.reflections, "Reflections have no bbox"

        # Select only those reflections used in refinement
        selection = self.reflections.get_flags(self.reflections.flags.reference_spot)
        if selection.count(True) == 0:
            raise RuntimeError("No reference reflections given")
        self.reflections = self.reflections.select(selection)

        # Compute the block size and jobs
        self.compute_blocks()
        self.compute_jobs()
        self.reflections = split_partials_over_boundaries(
            self.reflections, self.params.integration.block.size
        )

        # Create the reflection manager
        self.manager = SimpleReflectionManager(
            self.blocks, self.reflections, self.params.integration.mp.njobs
        )

    def task(self, index):
        """
        Get a task.
        """
        frames = self.manager.job(index)
        experiments = self.experiments
        reflections = self.manager.split(index)
        if len(reflections) == 0:
            logger.warning("No reflections in job %d ***", index)
            task = NullTask(index=index, reflections=reflections)
        else:
            task = ReferenceCalculatorJob(
                index=index,
                job=frames,
                experiments=experiments,
                reflections=reflections,
                params=self.params,
            )
        return task

    def tasks(self):
        """
        Iterate through the tasks.
        """
        for i in range(len(self)):
            yield self.task(i)

    def accumulate(self, result):
        """Accumulate the results."""
        self.manager.accumulate(result.index, result.reflections)

        if self.reference is None:
            self.reference = result.data
        else:
            self.reference.accumulate(result.data)

    def finalize(self):
        """
        Finalize the processing and finish.
        """
        # Check manager is finished
        assert self.manager.finished(), "Manager is not finished"

        # Set the reference profiles
        self.reference = self.reference.reference_profiles()

        self.finalized = True

    def result(self):
        """
        Return the result.

        :return: The result
        """
        assert self.finalized, "Manager is not finalized"
        return self.manager.data()

    def finished(self):
        """
        Return if all tasks have finished.

        :return: True/False all tasks have finished
        """
        return self.finalized and self.manager.finished()

    def __len__(self):
        """
        Return the number of tasks.

        :return: the number of tasks
        """
        return len(self.manager)

    def compute_max_block_size(self):
        """
        Compute the required memory
        """
        total_memory = psutil.virtual_memory().available
        max_memory_usage = self.params.integration.block.max_memory_usage
        assert max_memory_usage > 0.0, "maximum memory usage must be > 0"
        assert max_memory_usage <= 1.0, "maximum memory usage must be <= 1"
        limit_memory = int(math.floor(total_memory * max_memory_usage))
        return MultiThreadedReferenceProfiler.compute_max_block_size(
            self.experiments[0].imageset, max_memory_usage=limit_memory
        )

    def compute_blocks(self):
        """
        Compute the processing block size.
        """
        block = self.params.integration.block
        max_block_size = self.compute_max_block_size()
        if block.size in [Auto, "auto", "Auto"]:
            assert block.threshold > 0, "Threshold must be > 0"
            assert block.threshold <= 1.0, "Threshold must be < 1"
            nframes = sorted([b[5] - b[4] for b in self.reflections["bbox"]])
            cutoff = int(block.threshold * len(nframes))
            block_size = nframes[cutoff]
            if block_size > max_block_size:
                logger.warning(
                    "Computed block size (%s) > maximum block size (%s).",
                    block_size,
                    max_block_size,
                )
                logger.warning(
                    "Setting block size to maximum; some reflections may be partial"
                )
                block_size = max_block_size
        else:
            scan = self.experiments[0].scan
            if block.units == "radians":
                phi0, dphi = scan.get_oscillation(deg=False)
                block_size = int(math.ceil(block.size / dphi))
            elif block.units == "degrees":
                phi0, dphi = scan.get_oscillation()
                block_size = int(math.ceil(block.size / dphi))
            elif block.units == "frames":
                block_size = int(math.ceil(block.size))
            else:
                raise RuntimeError(f"Unknown block_size unit {block.units!r}")
            if block_size > max_block_size:
                raise RuntimeError(
                    """
          The requested block size (%s) is larger than the maximum allowable block
          size (%s). Either decrease the requested block size or increase the
          amount of available memory.
        """
                    % (block_size, max_block_size)
                )
        block.size = block_size
        block.units = "frames"

    def compute_jobs(self):
        imageset = self.experiments[0].imageset
        array_range = imageset.get_array_range()
        block = self.params.integration.block
        assert block.units == "frames"
        assert block.size > 0
        self.blocks = SimpleBlockList(array_range, block.size)
        assert len(self.blocks) > 0, "Invalid number of jobs"

    def summary(self):
        """
        Get a summary of the processing
        """
        # Compute the task table
        if self.experiments.all_stills():
            rows = [["#", "Group", "Frame From", "Frame To", "# Reflections"]]
            for i in range(len(self)):
                job = self.manager.job(i)
                group = job.index()
                f0, f1 = job.frames()
                n = self.manager.num_reflections(i)
                rows.append([str(i), str(group), str(f0), str(f1), str(n)])
        elif self.experiments.all_sequences():
            rows = [
                [
                    "#",
                    "Frame From",
                    "Frame To",
                    "Angle From",
                    "Angle To",
                    "# Reflections",
                ]
            ]
            for i in range(len(self)):
                f0, f1 = self.manager.job(i)
                scan = self.experiments[0].scan
                p0 = scan.get_angle_from_array_index(f0)
                p1 = scan.get_angle_from_array_index(f1)
                n = self.manager.num_reflections(i)
                rows.append([str(i), str(f0), str(f1), str(p0), str(p1), str(n)])
        else:
            raise RuntimeError("Experiments must be all sequences or all stills")

        # The job table
        task_table = tabulate(rows, headers="firstrow")

        # The format string
        if self.params.integration.block.size is None:
            block_size = "auto"
        else:
            block_size = str(self.params.integration.block.size)
        fmt = (
            "Processing reflections in the following blocks of images:\n"
            "\n"
            " block_size: %s %s\n"
            "\n"
            "%s\n"
        )
        return fmt % (block_size, self.params.integration.block.units, task_table)


def compute_required_memory(imageset, block_size):
    """
    Compute the required memory

    """
    return MultiThreadedIntegrator.compute_required_memory(imageset, block_size)


class ReferenceCalculatorProcessor:
    def __init__(self, experiments, reflections, params=None):
        from dials.util import pprint

        # Create the reference manager
        reference_manager = ReferenceCalculatorManager(experiments, reflections, params)

        # Print some output
        logger.info(reference_manager.summary())

        # Execute each task
        if params.integration.mp.njobs > 1:

            if params.integration.mp.method == "multiprocessing":
                _assert_enough_memory(
                    params.integration.mp.njobs
                    * compute_required_memory(
                        experiments[0].imageset, params.integration.block.size
                    ),
                    params.integration.block.max_memory_usage,
                )

            def process_output(result):
                for message in result[1]:
                    logger.log(message.levelno, message.msg)
                reference_manager.accumulate(result[0])

            multi_node_parallel_map(
                func=execute_parallel_task,
                iterable=reference_manager.tasks(),
                nproc=params.integration.mp.nproc,
                njobs=params.integration.mp.njobs,
                callback=process_output,
                cluster_method=params.integration.mp.method,
                preserve_order=True,
            )
        else:
            for task in reference_manager.tasks():
                result = task()
                reference_manager.accumulate(result)

        # Finalize the processing
        reference_manager.finalize()

        # Set the reflections and profiles
        self._reflections = reference_manager.result()
        self._profiles = reference_manager.reference

        # Write the profiles to file
        if params.integration.debug.reference.output:
            with open(params.integration.debug.reference.filename, "wb") as outfile:
                pickle.dump(self._profiles, outfile)

        # Print the profiles to the debug log
        for i in range(len(self._profiles)):
            logger.debug("")
            logger.debug("Reference Profiles for experiment %d", i)
            logger.debug("")
            reference = self._profiles[i].reference()
            for j in range(len(reference)):
                data = reference.data(j)
                logger.debug("Profile %d", j)
                if len(data) > 0:
                    logger.debug(pprint.profile3d(data))
                else:
                    logger.debug("** NO PROFILE **")

    def reflections(self):
        return self._reflections

    def profiles(self):
        return self._profiles


class IntegratorProcessor:
    def __init__(self, experiments, reflections, reference=None, params=None):

        # Create the reference manager
        integration_manager = IntegrationManager(
            experiments, reflections, reference, params
        )

        # Print some output
        logger.info(integration_manager.summary())

        # Execute each task
        if params.integration.mp.njobs > 1:

            if params.integration.mp.method == "multiprocessing":
                _assert_enough_memory(
                    params.integration.mp.njobs
                    * compute_required_memory(
                        experiments[0].imageset, params.integration.block.size
                    ),
                    params.integration.block.max_memory_usage,
                )

            def process_output(result):
                for message in result[1]:
                    logger.log(message.levelno, message.msg)
                integration_manager.accumulate(result[0])

            multi_node_parallel_map(
                func=execute_parallel_task,
                iterable=integration_manager.tasks(),
                nproc=params.integration.mp.nproc,
                njobs=params.integration.mp.njobs,
                callback=process_output,
                cluster_method=params.integration.mp.method,
                preserve_order=True,
            )
        else:
            for task in integration_manager.tasks():
                result = task()
                integration_manager.accumulate(result)

        # Finalize the processing
        integration_manager.finalize()

        # Set the reflections and profiles
        self._reflections = integration_manager.result()

    def reflections(self):
        return self._reflections


def split_partials_over_boundaries(reflections, block_size):
    """
    Split the reflections into partials or over job boundaries
    """

    # Get the block size and num frames
    _, _, _, _, z0, z1 = reflections["bbox"].parts()
    n_frames_of_bboxes = z1 - z0
    num_full = len(reflections)

    # See if any reflections need to be split
    refl_to_split_sel = n_frames_of_bboxes > block_size
    if refl_to_split_sel.count(True) > 0:

        # Get the subset of reflections to be split
        subset = reflections.select(refl_to_split_sel)
        newset = flex.reflection_table()
        for i in range(len(subset)):
            # get a copy of the reflection using selection
            item = subset.select(flex.size_t([i]))
            bbox = item["bbox"][0]
            size = bbox[5] - bbox[4]
            assert size > block_size
            nsplits = int(math.ceil(float(size) / float(block_size)))
            partsize = int(math.ceil(float(size) / float(nsplits)))
            bbox0 = bbox[0:4] + (bbox[4], bbox[4] + partsize)
            subset["bbox"][i] = bbox0  # set the updated bbox in the subset
            for _ in range(1, nsplits):
                # make z0 start at end of previous bbox0
                bbox0 = bbox0[0:4] + (bbox0[5], min(bbox0[5] + partsize, bbox[5]))
                assert bbox0[4] < bbox0[5]
                assert bbox0[5] <= bbox[5]
                # set the split bbox as the updated bbox in the reflection copy
                item["bbox"][0] = bbox0
                # and add to newset as separate reflection
                newset.extend(item)

        # Update the subset of trimmed reflections and add the new set of partials
        reflections.set_selected(refl_to_split_sel, subset)
        reflections.extend(newset)

    # Print some info
    num_after_splitting = len(reflections)
    assert num_after_splitting >= num_full, "Invalid number of partials"
    if num_after_splitting > num_full:
        num_split = num_after_splitting - num_full
        logger.info(" Split %d reflections\n", num_split)

    return reflections
