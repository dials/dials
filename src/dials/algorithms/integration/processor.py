from __future__ import annotations

import itertools
import logging
import math
import os
import pathlib
from time import time

import psutil

import boost_adaptbx.boost.python
import libtbx

import dials.algorithms.integration
import dials.util
import dials.util.log
from dials.array_family import flex
from dials.model.data import make_image
from dials.util import tabulate
from dials.util.log import rehandle_cached_records
from dials.util.mp import available_cores, multi_node_parallel_map
from dials_algorithms_integration_integrator_ext import (
    Executor,
    Group,
    GroupList,
    Job,
    JobList,
    ReflectionManager,
    ReflectionManagerPerImage,
    ShoeboxProcessor,
)

try:
    import resource
except ImportError:
    # resource does not exist on non-Linux, so can't float the import
    resource = None

__all__ = [
    "Block",
    "build_processor",
    "Debug",
    "Executor",
    "Group",
    "GroupList",
    "Job",
    "job",
    "JobList",
    "Lookup",
    "MultiProcessing",
    "NullTask",
    "Parameters",
    "Processor2D",
    "Processor3D",
    "ProcessorFlat3D",
    "ProcessorSingle2D",
    "ProcessorStills",
    "ReflectionManager",
    "ReflectionManagerPerImage",
    "Shoebox",
    "ShoeboxProcessor",
    "Task",
]

logger = logging.getLogger(__name__)


def assess_available_memory(params):

    # Obtain information about system memory
    available_memory = psutil.virtual_memory().available
    available_swap = psutil.swap_memory().free

    # https://htcondor.readthedocs.io/en/latest/users-manual/services-for-jobs.html#extra-environment-variables-htcondor-sets-for-jobs
    condor_job_ad = os.environ.get("_CONDOR_JOB_AD")
    if condor_job_ad:
        try:
            job_classad = dials.util.parse_htcondor_job_classad(
                pathlib.Path(condor_job_ad)
            )
        except Exception as e:
            logger.error(
                f"Error parsing _CONDOR_JOB_AD {condor_job_ad}: {e}",
                exc_info=True,
            )
        else:
            if job_classad.memory_provisioned:
                # Convert MB to bytes
                available_memory = min(
                    available_memory, job_classad.memory_provisioned * 1024**2
                )

    available_incl_swap = available_memory + available_swap
    available_limit = available_incl_swap * params.block.max_memory_usage
    available_immediate_limit = available_memory * params.block.max_memory_usage

    # Compile a memory report
    report = [
        "Memory situation report:",
    ]

    def _report(description, numbytes):
        report.append(f"  {description:<50}:{numbytes/1e9:5.1f} GB")

    _report("Available system memory (excluding swap)", available_memory)
    _report("Available swap memory", available_swap)
    _report("Available system memory (including swap)", available_incl_swap)
    _report("Maximum memory for processing (including swap)", available_limit)
    _report(
        "Maximum memory for processing (excluding swap)",
        available_immediate_limit,
    )

    # Check if a ulimit applies
    # Note that resource may be None on non-Linux platforms.
    # We can't use psutil as platform-independent solution in this instance due to
    # https://github.com/conda-forge/psutil-feedstock/issues/47
    rlimit = getattr(resource, "RLIMIT_VMEM", getattr(resource, "RLIMIT_AS", None))
    if rlimit:
        try:
            ulimit = resource.getrlimit(rlimit)[0]
            if ulimit <= 0 or ulimit > (2**62):
                report.append("  no memory ulimit set")
            else:
                ulimit_used = psutil.Process().memory_info().rss
                _report("Memory ulimit detected", ulimit)
                _report("Memory ulimit in use", ulimit_used)
                available_memory = max(0, min(available_memory, ulimit - ulimit_used))
                available_incl_swap = max(
                    0, min(available_incl_swap, ulimit - ulimit_used)
                )
                available_immediate_limit = (
                    available_memory * params.block.max_memory_usage
                )
                _report("Available system memory (limited)", available_memory)
                _report(
                    "Available system memory (incl. swap; limited)",
                    available_incl_swap,
                )
                _report(
                    "Maximum memory for processing (exc. swap; limited)",
                    available_immediate_limit,
                )
        except Exception as e:
            logger.debug(
                "Could not obtain ulimit values due to %s", str(e), exc_info=True
            )

    return available_immediate_limit, available_incl_swap, report


def _average_bbox_size(reflections):
    """Calculate the average bbox size for debugging"""

    bbox = reflections["bbox"]
    sel = flex.random_selection(len(bbox), min(len(bbox), 1000))
    subset_bbox = bbox.select(sel)
    xmin, xmax, ymin, ymax, zmin, zmax = subset_bbox.parts()
    xsize = flex.mean((xmax - xmin).as_double())
    ysize = flex.mean((ymax - ymin).as_double())
    zsize = flex.mean((zmax - zmin).as_double())
    return xsize, ysize, zsize


@boost_adaptbx.boost.python.inject_into(Executor)
class _:
    @staticmethod
    def __getinitargs__():
        return ()


class _Job:
    def __init__(self):
        self.index = 0
        self.nthreads = 1


job = _Job()


class MultiProcessing:
    """
    Multi processing parameters
    """

    def __init__(self):
        self.method = "multiprocessing"
        self.nproc = 1
        self.njobs = 1
        self.nthreads = 1
        self.n_subset_split = None

    def update(self, other):
        self.method = other.method
        self.nproc = other.nproc
        self.njobs = other.njobs
        self.nthreads = other.nthreads
        self.n_subset_split = other.n_subset_split


class Lookup:
    """
    Lookup parameters
    """

    def __init__(self):
        self.mask = None

    def update(self, other):
        self.mask = other.mask


class Block:
    """
    Block parameters
    """

    def __init__(self):
        self.size = libtbx.Auto
        self.units = "degrees"
        self.threshold = 0.99
        self.force = False
        self.max_memory_usage = 0.90

    def update(self, other):
        self.size = other.size
        self.units = other.units
        self.threshold = other.threshold
        self.force = other.force
        self.max_memory_usage = other.max_memory_usage


class Shoebox:
    """
    Shoebox parameters
    """

    def __init__(self):
        self.flatten = False
        self.partials = False

    def update(self, other):
        self.flatten = other.flatten
        self.partials = other.partials


class Debug:
    """
    Debug parameters
    """

    def __init__(self):
        self.output = False
        self.select = None
        self.split_experiments = True
        self.separate_files = True

    def update(self, other):
        self.output = other.output
        self.select = other.select
        self.split_experiments = other.split_experiments
        self.separate_files = other.separate_files


class Parameters:
    """
    Class to handle parameters for the processor
    """

    def __init__(self):
        """
        Initialize the parameters
        """
        self.mp = MultiProcessing()
        self.lookup = Lookup()
        self.block = Block()
        self.shoebox = Shoebox()
        self.debug = Debug()

    def update(self, other):
        """
        Update the parameters
        """
        self.mp.update(other.mp)
        self.lookup.update(other.lookup)
        self.block.update(other.block)
        self.shoebox.update(other.shoebox)
        self.debug.update(other.debug)


def execute_parallel_task(task):
    """
    Helper function to run things on cluster
    """

    dials.util.log.config_simple_cached()
    result = task()
    handlers = logging.getLogger("dials").handlers
    assert len(handlers) == 1, "Invalid number of logging handlers"
    return result, handlers[0].records


class _Processor:
    """Processor interface class."""

    def __init__(self, manager):
        """
        Initialise the processor.

        The processor requires a manager class implementing the _Manager interface.
        This class executes all the workers in separate threads and accumulates the
        results to expose to the user.

        :param manager: The processing manager
        :param params: The phil parameters
        """
        self.manager = manager

    @property
    def executor(self):
        """
        Get the executor

        :return: The executor
        """
        return self.manager.executor

    @executor.setter
    def executor(self, function):
        """
        Set the executor

        :param function: The executor
        """
        self.manager.executor = function

    def process(self):
        """
        Do all the processing tasks.

        :return: The processing results
        """
        start_time = time()
        self.manager.initialize()
        mp_method = self.manager.params.mp.method
        mp_njobs = self.manager.params.mp.njobs
        mp_nproc = self.manager.params.mp.nproc

        assert mp_nproc > 0, "Invalid number of processors"
        if mp_nproc * mp_njobs > len(self.manager):
            mp_nproc = min(mp_nproc, len(self.manager))
            mp_njobs = int(math.ceil(len(self.manager) / mp_nproc))
        logger.info(self.manager.summary())
        if mp_njobs > 1:
            assert mp_method != "none" and mp_method is not None
            logger.info(
                " Using %s with %d parallel job(s) and %d processes per node\n",
                mp_method,
                mp_njobs,
                mp_nproc,
            )
        else:
            logger.info(" Using multiprocessing with %d parallel job(s)\n", mp_nproc)

        if mp_njobs * mp_nproc > 1:

            def process_output(result):
                rehandle_cached_records(result[1])
                self.manager.accumulate(result[0])

            multi_node_parallel_map(
                func=execute_parallel_task,
                iterable=list(self.manager.tasks()),
                njobs=mp_njobs,
                nproc=mp_nproc,
                callback=process_output,
                cluster_method=mp_method,
                preserve_order=True,
            )
        else:
            for task in self.manager.tasks():
                self.manager.accumulate(task())
        self.manager.finalize()
        end_time = time()
        self.manager.time.user_time = end_time - start_time
        result1, result2 = self.manager.result()
        return result1, result2, self.manager.time


class _ProcessorRot(_Processor):
    """Processor interface class for rotation data only."""

    def __init__(self, experiments, manager):
        """
        Initialise the processor.

        The processor requires a manager class implementing the _Manager interface.
        This class executes all the workers in separate threads and accumulates the
        results to expose to the user.

        :param manager: The processing manager
        """
        # Ensure we have the correct type of data
        if not experiments.all_sequences():
            raise RuntimeError(
                """
        An inappropriate processing algorithm may have been selected!
         Trying to perform rotation processing when not all experiments
         are indicated as rotation experiments.
      """
            )

        super().__init__(manager)


class NullTask:
    """
    A class to perform a null task.
    """

    def __init__(self, index, reflections):
        """
        Initialise the task

        :param index: The index of the processing job
        :param experiments: The list of experiments
        :param reflections: The list of reflections
        """
        self.index = index
        self.reflections = reflections

    def __call__(self):
        """
        Do the processing.

        :return: The processed data
        """
        return dials.algorithms.integration.Result(
            index=self.index,
            reflections=self.reflections,
            data=None,
            read_time=0,
            extract_time=0,
            process_time=0,
            total_time=0,
        )


class Task:
    """
    A class to perform a processing task.
    """

    def __init__(self, index, job, experiments, reflections, params, executor=None):
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
        assert executor is not None, "No executor given"
        assert len(reflections) > 0, "Zero reflections given"
        self.index = index
        self.job = job
        self.experiments = experiments
        self.reflections = reflections
        self.params = params
        self.executor = executor

    def __call__(self):
        """
        Do the processing.

        :return: The processed data
        """
        # Get the start time
        start_time = time()

        # Set the global process ID
        job.index = self.index

        # Check all reflections have same imageset and get it
        exp_id = list(set(self.reflections["id"]))
        imageset = self.experiments[exp_id[0]].imageset
        for i in exp_id[1:]:
            assert (
                self.experiments[i].imageset == imageset
            ), "Task can only handle 1 imageset"

        # Get the sub imageset
        frame0, frame1 = self.job

        try:
            allowed_range = imageset.get_array_range()
        except Exception:
            allowed_range = 0, len(imageset)

        try:
            # range increasing
            assert frame0 < frame1

            # within an increasing range
            assert allowed_range[1] > allowed_range[0]

            # we are processing data which is within range
            assert frame0 >= allowed_range[0]
            assert frame1 <= allowed_range[1]

            # I am 99% sure this is implied by all the code above
            assert (frame1 - frame0) <= len(imageset)
            if len(imageset) > 1:
                imageset = imageset[frame0:frame1]
        except Exception as e:
            raise RuntimeError(f"Programmer Error: bad array range: {e}")

        try:
            frame0, frame1 = imageset.get_array_range()
        except Exception:
            frame0, frame1 = (0, len(imageset))

        self.executor.initialize(frame0, frame1, self.reflections)

        # Set the shoeboxes (don't allocate)
        self.reflections["shoebox"] = flex.shoebox(
            self.reflections["panel"],
            self.reflections["bbox"],
            allocate=False,
            flatten=self.params.shoebox.flatten,
        )

        # Create the processor
        processor = ShoeboxProcessor(
            self.reflections,
            len(imageset.get_detector()),
            frame0,
            frame1,
            self.params.debug.output,
        )

        # Loop through the imageset, extract pixels and process reflections
        read_time = 0.0
        for i in range(len(imageset)):
            st = time()
            image = imageset.get_corrected_data(i)
            if imageset.is_marked_for_rejection(i):
                mask = tuple(flex.bool(im.accessor(), False) for im in image)
            else:
                mask = imageset.get_mask(i)
                if self.params.lookup.mask is not None:
                    assert len(mask) == len(
                        self.params.lookup.mask
                    ), "Mask/Image are incorrect size %d %d" % (
                        len(mask),
                        len(self.params.lookup.mask),
                    )
                    mask = tuple(
                        m1 & m2 for m1, m2 in zip(self.params.lookup.mask, mask)
                    )

            read_time += time() - st
            processor.next(make_image(image, mask), self.executor)
            del image
            del mask
        assert processor.finished(), "Data processor is not finished"

        # Optionally save the shoeboxes
        if self.params.debug.output and self.params.debug.separate_files:
            output = self.reflections
            if self.params.debug.select is not None:
                output = output.select(self.params.debug.select(output))
            if self.params.debug.split_experiments:
                output = output.split_by_experiment_id()
                for table in output:
                    i = table["id"][0]
                    table.as_file("shoeboxes_%d_%d.refl" % (self.index, i))
            else:
                output.as_file("shoeboxes_%d.refl" % self.index)

        # Delete the shoeboxes
        if self.params.debug.separate_files or not self.params.debug.output:
            del self.reflections["shoebox"]

        # Finalize the executor
        self.executor.finalize()

        # Return the result
        return dials.algorithms.integration.Result(
            index=self.index,
            reflections=self.reflections,
            data=self.executor.data(),
            read_time=read_time,
            extract_time=processor.extract_time(),
            process_time=processor.process_time(),
            total_time=time() - start_time,
        )


class _Manager:
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

        # Initialise the callbacks
        self.executor = None

        # Save some data
        self.experiments = experiments
        self.reflections = reflections

        # Other data
        self.data = {}

        # Save some parameters
        self.params = params

        # Set the finalized flag to False
        self.finalized = False

        # Initialise the timing information
        self.time = dials.algorithms.integration.TimingInfo()

    def initialize(self):
        """
        Initialise the processing
        """
        # Get the start time
        start_time = time()

        # Ensure the reflections contain bounding boxes
        assert "bbox" in self.reflections, "Reflections have no bbox"

        if self.params.mp.nproc is libtbx.Auto:
            self.params.mp.nproc = available_cores()
            logger.info(f"Setting nproc={self.params.mp.nproc}")

        # Compute the block size and processors
        self.compute_jobs()
        self.split_reflections()
        self.compute_processors()

        # Create the reflection manager
        self.manager = ReflectionManager(self.jobs, self.reflections)

        # Set the initialization time
        self.time.initialize = time() - start_time

    def task(self, index):
        """
        Get a task.
        """
        job = self.manager.job(index)
        frames = job.frames()
        expr_id = job.expr()
        assert expr_id[1] > expr_id[0], "Invalid experiment id"
        assert expr_id[0] >= 0, "Invalid experiment id"
        assert expr_id[1] <= len(self.experiments), "Invalid experiment id"
        experiments = self.experiments  # [expr_id[0]:expr_id[1]]
        reflections = self.manager.split(index)
        if len(reflections) == 0:
            logger.warning("No reflections in job %d ***", index)
            task = NullTask(index=index, reflections=reflections)
        else:
            task = Task(
                index=index,
                job=frames,
                experiments=experiments,
                reflections=reflections,
                params=self.params,
                executor=self.executor,
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
        self.data[result.index] = result.data
        self.manager.accumulate(result.index, result.reflections)
        self.time.read += result.read_time
        self.time.extract += result.extract_time
        self.time.process += result.process_time
        self.time.total += result.total_time

    def finalize(self):
        """
        Finalize the processing and finish.
        """
        # Get the start time
        start_time = time()

        # Check manager is finished
        assert self.manager.finished(), "Manager is not finished"

        # Update the time and finalized flag
        self.time.finalize = time() - start_time
        self.finalized = True

    def result(self):
        """
        Return the result.

        :return: The result
        """
        assert self.finalized, "Manager is not finalized"
        return self.manager.data(), self.data

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

    def compute_jobs(self):
        """
        Sets up a JobList() object in self.jobs
        """

        if self.params.block.size == libtbx.Auto:
            if (
                self.params.mp.nproc * self.params.mp.njobs == 1
                and not self.params.debug.output
                and not self.params.block.force
            ):
                self.params.block.size = None

        # calculate the block overlap based on the size of bboxes in the data
        # calculate once here rather than repeated in the loop below
        block_overlap = 0
        if self.params.block.size is not None:
            assert self.params.block.threshold > 0, "Threshold must be > 0"
            assert self.params.block.threshold <= 1.0, "Threshold must be < 1"
            frames_per_refl = sorted([b[5] - b[4] for b in self.reflections["bbox"]])
            cutoff = int(self.params.block.threshold * len(frames_per_refl))
            block_overlap = frames_per_refl[cutoff]

        groups = itertools.groupby(
            range(len(self.experiments)),
            lambda x: (id(self.experiments[x].imageset), id(self.experiments[x].scan)),
        )
        self.jobs = JobList()
        for key, indices in groups:
            indices = list(indices)
            i0 = indices[0]
            i1 = indices[-1] + 1
            expr = self.experiments[i0]
            scan = expr.scan
            imgs = expr.imageset
            array_range = (0, len(imgs))
            if scan is not None:
                assert len(imgs) >= len(scan), "Invalid scan range"
                array_range = scan.get_array_range()

            if self.params.block.size is None:
                block_size_frames = array_range[1] - array_range[0]
            elif self.params.block.size == libtbx.Auto:
                # auto determine based on nframes and overlap
                nframes = array_range[1] - array_range[0]
                nblocks = self.params.mp.nproc * self.params.mp.njobs
                # want data to be split into n blocks with overlaps
                # i.e. [x, overlap, y, overlap, y, overlap, ....,y,  overlap, x]
                # blocks are x + overlap, or overlap + y + overlap.
                x = (nframes - block_overlap) / nblocks
                block_size = int(math.ceil(x + block_overlap))
                # increase the block size to be at least twice the overlap, in
                # case the overlap is large e.g. if high mosaicity.
                block_size_frames = max(block_size, 2 * block_overlap)
            elif self.params.block.units == "radians":
                _, dphi = scan.get_oscillation(deg=False)
                block_size_frames = int(math.ceil(self.params.block.size / dphi))
                # if the specified block size is lower than the overlap,
                # reduce the overlap to be half of the block size.
                block_overlap = min(block_overlap, int(block_size_frames // 2))
            elif self.params.block.units == "degrees":
                _, dphi = scan.get_oscillation()
                block_size_frames = int(math.ceil(self.params.block.size / dphi))
                # if the specified block size is lower than the overlap,
                # reduce the overlap to be half of the block size.
                block_overlap = min(block_overlap, int(block_size_frames // 2))
            elif self.params.block.units == "frames":
                block_size_frames = int(math.ceil(self.params.block.size))
                block_overlap = min(block_overlap, int(block_size_frames // 2))
            else:
                raise RuntimeError(
                    f"Unknown block_size units {self.params.block.units!r}"
                )
            self.jobs.add(
                (i0, i1),
                array_range,
                block_size_frames,
                block_overlap,
            )
        assert len(self.jobs) > 0, "Invalid number of jobs"

    def split_reflections(self):
        """
        Split the reflections into partials or over job boundaries
        """

        # Optionally split the reflection table into partials, otherwise,
        # split over job boundaries
        if self.params.shoebox.partials:
            num_full = len(self.reflections)
            self.reflections.split_partials()
            num_partial = len(self.reflections)
            assert num_partial >= num_full, "Invalid number of partials"
            if num_partial > num_full:
                logger.info(
                    " Split %d reflections into %d partial reflections\n",
                    num_full,
                    num_partial,
                )
        else:
            num_full = len(self.reflections)
            self.jobs.split(self.reflections)
            num_partial = len(self.reflections)
            assert num_partial >= num_full, "Invalid number of partials"
            if num_partial > num_full:
                num_split = num_partial - num_full
                logger.info(
                    " Split %d reflections overlapping job boundaries\n", num_split
                )

        # Compute the partiality
        self.reflections.compute_partiality(self.experiments)

    def compute_processors(self):
        """
        Compute the number of processors
        """

        # Get the maximum shoebox memory to estimate memory use for one process
        memory_required_per_process = flex.max(
            self.jobs.shoebox_memory(self.reflections, self.params.shoebox.flatten)
        )

        (
            available_immediate_limit,
            available_incl_swap,
            report,
        ) = assess_available_memory(self.params)

        report.append(
            f"  {'Memory required per process':50}:{memory_required_per_process/1e9:5.1f} GB"
        )

        output_level = logging.INFO

        # Limit the number of parallel processes by amount of available memory
        if self.params.mp.method == "multiprocessing" and self.params.mp.nproc > 1:

            # Compute expected memory usage and warn if not enough
            njobs = available_immediate_limit / memory_required_per_process
            if njobs >= self.params.mp.nproc:
                # There is enough memory. Take no action
                pass
            elif njobs >= 1:
                # There is enough memory to run, but not as many processes as requested
                output_level = logging.WARNING
                report.append(
                    "Reducing number of processes from %d to %d due to memory constraints."
                    % (self.params.mp.nproc, int(njobs))
                )
                self.params.mp.nproc = int(njobs)
            elif (
                available_incl_swap * self.params.block.max_memory_usage
                >= memory_required_per_process
            ):
                # There is enough memory to run, but only if we count swap.
                output_level = logging.WARNING
                report.append(
                    "Reducing number of processes from %d to 1 due to memory constraints."
                    % self.params.mp.nproc,
                )
                report.append("Running this process will rely on swap usage!")
                self.params.mp.nproc = 1
            else:
                # There is not enough memory to run
                output_level = logging.ERROR

        report.append("")
        logger.log(output_level, "\n".join(report))

        if output_level >= logging.ERROR:
            raise MemoryError(
                """
          Not enough memory to run integration jobs.  This could be caused by a
          highly mosaic crystal model.  Possible solutions include increasing the
          percentage of memory allowed for shoeboxes or decreasing the block size.
          The average shoebox size is %d x %d pixels x %d images - is your crystal
          really this mosaic?
          """
                % _average_bbox_size(self.reflections)
            )

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
                    "Group",
                    "Frame From",
                    "Frame To",
                    "Angle From",
                    "Angle To",
                    "# Reflections",
                ]
            ]
            for i in range(len(self)):
                job = self.manager.job(i)
                group = job.index()
                expr = job.expr()
                f0, f1 = job.frames()
                scan = self.experiments[expr[0]].scan
                p0 = scan.get_angle_from_array_index(f0)
                p1 = scan.get_angle_from_array_index(f1)
                n = self.manager.num_reflections(i)
                rows.append(
                    [str(i), str(group), str(f0 + 1), str(f1), str(p0), str(p1), str(n)]
                )
        else:
            raise RuntimeError("Experiments must be all sequences or all stills")

        # The job table
        task_table = tabulate(rows, headers="firstrow")

        # The format string
        if self.params.block.size is None:
            block_size = "auto"
        else:
            block_size = str(self.params.block.size)
        return (
            "Processing reflections in the following blocks of images:\n\n"
            " block_size: {} {}\n\n{}\n"
        ).format(
            block_size,
            "" if block_size in ("auto", "Auto") else self.params.block.units,
            task_table,
        )


class Processor3D(_ProcessorRot):
    """Top level processor for 3D processing."""

    def __init__(self, experiments, reflections, params):
        """Initialise the manager and the processor."""

        # Set some parameters
        params.shoebox.partials = False
        params.shoebox.flatten = False

        # Create the processing manager
        manager = _Manager(experiments, reflections, params)

        # Initialise the processor
        super().__init__(experiments, manager)


class ProcessorFlat3D(_ProcessorRot):
    """Top level processor for flat 3D processing."""

    def __init__(self, experiments, reflections, params):
        """Initialise the manager and the processor."""

        # Set some parameters
        params.shoebox.flatten = True
        params.shoebox.partials = False

        # Create the processing manager
        manager = _Manager(experiments, reflections, params)

        # Initialise the processor
        super().__init__(experiments, manager)


class Processor2D(_ProcessorRot):
    """Top level processor for 2D processing."""

    def __init__(self, experiments, reflections, params):
        """Initialise the manager and the processor."""

        # Set some parameters
        params.shoebox.partials = True

        # Create the processing manager
        manager = _Manager(experiments, reflections, params)

        # Initialise the processor
        super().__init__(experiments, manager)


class ProcessorSingle2D(_ProcessorRot):
    """Top level processor for still image processing."""

    def __init__(self, experiments, reflections, params):
        """Initialise the manager and the processor."""

        # Set some of the parameters
        params.block.size = 1
        params.block.units = "frames"
        params.shoebox.partials = True
        params.shoebox.flatten = False

        # Create the processing manager
        manager = _Manager(experiments, reflections, params)

        # Initialise the processor
        super().__init__(experiments, manager)


class ProcessorStills(_Processor):
    """Top level processor for still image processing."""

    def __init__(self, experiments, reflections, params):
        """Initialise the manager and the processor."""

        # Set some parameters
        params.block.size = 1
        params.block.units = "frames"
        params.shoebox.partials = False
        params.shoebox.flatten = False

        # Ensure we have the correct type of data
        if not experiments.all_stills():
            raise RuntimeError(
                """
        An inappropriate processing algorithm may have been selected!
         Trying to perform stills processing when not all experiments
         are indicated as stills experiments.
      """
            )

        # Create the processing manager
        manager = _Manager(experiments, reflections, params)

        # Initialise the processor
        super().__init__(manager)


def build_processor(Class, experiments, reflections, params=None):
    """
    A function to simplify building the processor

    :param Class: The input class
    :param experiments: The input experiments
    :param reflections: The reflections
    :param params: Optional input parameters
    """
    _params = Parameters()
    if params is not None:
        _params.update(params)

    return Class(experiments, reflections, _params)
