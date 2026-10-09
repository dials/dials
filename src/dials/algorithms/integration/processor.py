from __future__ import annotations

import bisect
import itertools
import logging
import math
import os
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
from dials.util.mp import multi_node_parallel_map
from dials.util.system import CPU_COUNT, MEMORY_LIMIT
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
    "OnePassProcessor3D",
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
            assert self.experiments[i].imageset == imageset, (
                "Task can only handle 1 imageset"
            )

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
                # Slice imageset as a 0-based array
                index0 = frame0 - allowed_range[0]
                index1 = frame1 - allowed_range[0]
                imageset = imageset[index0:index1]
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
                    assert len(mask) == len(self.params.lookup.mask), (
                        "Mask/Image are incorrect size %d %d"
                        % (
                            len(mask),
                            len(self.params.lookup.mask),
                        )
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
            self.params.mp.nproc = CPU_COUNT
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

        # Obtain information about system memory
        available_memory = MEMORY_LIMIT
        available_limit = available_memory
        if self.params.block.max_memory_usage is not None:
            available_limit *= self.params.block.max_memory_usage

        # Get the maximum shoebox memory to estimate memory use for one process
        required_shoebox_memory = self.required_shoebox_memory()

        # Get the current memory usage
        current_memory_usage = psutil.Process(os.getpid()).memory_info().rss

        memory_required_per_process = required_shoebox_memory + current_memory_usage

        # Compile a memory report
        report = ["Memory situation report:"]

        def _report(description, numbytes):
            report.append(f"  {description:<50}: {numbytes / 1e6:5.1f} MB")

        _report("Available system memory", available_memory)
        _report("Maximum memory for processing", available_limit)
        _report("Current memory usage", current_memory_usage)
        _report("Memory required for shoeboxes", required_shoebox_memory)
        _report("Memory required per process", memory_required_per_process)

        output_level = logging.INFO

        # Limit the number of parallel processes by amount of available memory
        if (
            self.params.mp.method == "multiprocessing"
            and self.params.mp.nproc > 1
            and self.params.block.max_memory_usage is not None
        ):
            # Compute expected memory usage and warn if not enough
            njobs = available_limit / memory_required_per_process
            if njobs >= self.params.mp.nproc:
                # There is enough memory. Take no action
                pass
            elif njobs >= 1:
                # There is enough memory to run, but not as many processes as requested
                output_level = logging.WARNING
                report.append(
                    f"Reducing number of processes from {self.params.mp.nproc} to "
                    f"{int(njobs)} due to memory constraints."
                )
                self.params.mp.nproc = int(njobs)
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

    def required_shoebox_memory(self):
        """
        The most memory any one job needs for shoeboxes
        """
        return flex.max(
            self.jobs.shoebox_memory(self.reflections, self.params.shoebox.flatten)
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


class _OnePassJob:
    """
    A job for one-pass integration: the reflections it owns, the reference
    spots it reads only to model profiles, and the frames they span.
    """

    def __init__(self, group, expr, frames, rows, learn_only, owned_cells):
        self.group = group
        self.expr = expr
        self.frames = frames
        self.rows = rows
        self.learn_only = learn_only
        self.owned_cells = owned_cells


class _OnePassManager(_Manager):
    """
    Book-keeping for integrating in one pass with several jobs.

    The reference profiles of each imageset are divided by position in the scan
    into contiguous groups, one for each job. A job owns the reflections fitted
    against its profiles, and also reads every reference spot that could be
    added to them, wherever it is owned, so that it can model them completely
    itself: no profile is shared between jobs. A reference profile is the sum of
    its contributions in the order the images are read, so each job's profiles
    are exactly those of a single job reading every image. Each job reads the
    frames its reflections span, so jobs overlap where they share reference
    spots, and no reflection is split between jobs.
    """

    LEARN_ONLY = "one_pass.learn_only"

    # A job reads the reference spots up to one block of profiles beyond those
    # it owns on each side, so owning fewer than this many blocks would mostly
    # add reading rather than parallelism
    MIN_BLOCKS_PER_JOB = 3

    def __init__(self, experiments, reflections, params, profile_modeller):
        super().__init__(experiments, reflections, params)
        self.profile_modeller = profile_modeller

    def initialize(self):
        """
        Initialise the processing
        """
        start_time = time()
        assert "bbox" in self.reflections, "Reflections have no bbox"
        if self.params.mp.nproc is libtbx.Auto:
            self.params.mp.nproc = CPU_COUNT
            logger.info(f"Setting nproc={self.params.mp.nproc}")
        self.compute_jobs()
        # As JobList.split sets it when no reflection is split
        self.reflections["partial_id"] = flex.size_t_range(len(self.reflections))
        self.reflections.compute_partiality(self.experiments)
        self.compute_processors()
        self.finished_jobs = [False] * len(self.one_pass_jobs)
        self.time.initialize = time() - start_time

    def compute_jobs(self):
        """
        Divide each imageset's reference profiles between the jobs and work out
        the reflections and frames of each.
        """
        reflections = self.reflections
        flags = reflections.flags
        processed = ~reflections.get_flags(flags.dont_integrate)
        reference = reflections.get_flags(flags.reference_spot) & processed
        ids = reflections["id"]
        z = reflections["xyzcal.px"].parts()[2]
        bbox = reflections["bbox"]
        nblocks = self.params.mp.nproc * self.params.mp.njobs

        groups = itertools.groupby(
            range(len(self.experiments)),
            lambda x: (id(self.experiments[x].imageset), id(self.experiments[x].scan)),
        )
        self.one_pass_jobs = []
        self.fitting_cell = flex.int(len(reflections), -1)
        for group, (_, indices) in enumerate(groups):
            indices = list(indices)
            i0, i1 = indices[0], indices[-1] + 1
            array_range = self.experiments[i0].scan.get_array_range()

            # The scan positions of the profiles, divided into contiguous groups
            positions = sorted(
                {
                    self.profile_modeller[i].coord(c)[2]
                    for i in indices
                    for c in range(len(self.profile_modeller[i]))
                }
            )
            n = max(1, min(nblocks, len(positions) // self.MIN_BLOCKS_PER_JOB))
            chunks = [
                positions[k * len(positions) // n : (k + 1) * len(positions) // n]
                for k in range(n)
            ]
            chunk_of = {p: k for k, chunk in enumerate(chunks) for p in chunk}
            boundaries = [(chunks[k][-1] + chunks[k + 1][0]) / 2 for k in range(n - 1)]

            # Each reflection belongs to the job owning the profile it is fitted
            # against; one with none, by its position in the scan
            owner = flex.int(len(reflections), -1)
            owned_cells = [{} for _ in range(n)]
            contributes = [flex.bool(len(reflections), False) for _ in range(n)]
            for i in indices:
                modeller = self.profile_modeller[i]
                chunk_of_cell = [
                    chunk_of[modeller.coord(c)[2]] for c in range(len(modeller))
                ]
                selection = (ids == i) & processed
                rows = selection.iselection()
                cells = modeller.fitting_cells(reflections.select(rows))
                self.fitting_cell.set_selected(rows, cells)
                for row, cell in zip(rows, cells):
                    if cell >= 0:
                        owner[row] = chunk_of_cell[cell]
                    else:
                        owner[row] = bisect.bisect(boundaries, z[row])
                reference_rows = (selection & reference).iselection()
                reference_table = reflections.select(reference_rows)
                for k in range(n):
                    mask = flex.bool([c == k for c in chunk_of_cell])
                    owned_cells[k][i] = mask.iselection()
                    contributes[k].set_selected(
                        reference_rows,
                        modeller.contributes_to(reference_table, mask),
                    )

            for k in range(n):
                owned = owner == k
                rows = (owned | contributes[k]).iselection()
                if len(rows) == 0:
                    continue
                learn_only = ~owned.select(rows)
                b = bbox.select(rows).parts()
                frames = (flex.min(b[4]), flex.max(b[5]))
                assert frames[0] >= array_range[0] and frames[1] <= array_range[1]
                self.one_pass_jobs.append(
                    _OnePassJob(
                        group, (i0, i1), frames, rows, learn_only, owned_cells[k]
                    )
                )
        assert len(self.one_pass_jobs) > 0, "Invalid number of jobs"

    def required_shoebox_memory(self):
        """
        The most memory any one job needs for shoeboxes, and for the transformed
        shoeboxes it holds until their reference profiles are finalized.
        """
        bbox = self.reflections["bbox"]
        flags = self.reflections.flags
        reference = self.reflections.get_flags(flags.reference_spot)
        required = 0
        for job in self.one_pass_jobs:
            frame0, frame1 = job.frames
            nframes = frame1 - frame0
            usage = [0] * (nframes + 1)
            x0, x1, y0, y1, z0, z1 = (
                a.as_numpy_array() for a in bbox.select(job.rows).parts()
            )

            # Shoeboxes, from their first frame to their last, as JobList counts
            nbytes = (x1 - x0) * (y1 - y0) * (z1 - z0) * 12
            for first, last, n in zip(z0 - frame0, z1 - frame0, nbytes):
                usage[first] += int(n)
                usage[last] -= int(n)

            # Transformed shoeboxes, from their last frame until the reference
            # profile they are fitted against is finalized
            table = self.reflections.select(job.rows)
            owned = ~job.learn_only
            for i in sorted(set(table["id"])):
                experiment = self.experiments[i]
                grid_size = experiment.profile.params.gaussian_rs.fitting.grid_size
                transform = (2 * grid_size + 1) ** 3 * (8 + 8 + 1)
                candidates = table.select(
                    (table["id"] == i) & reference.select(job.rows)
                )
                deadlines = self.profile_modeller[i].learning_deadlines(candidates)
                rows = ((table["id"] == i) & owned).iselection()
                cells = self.fitting_cell.select(job.rows.select(rows))
                for row, cell in zip(rows, cells):
                    if cell < 0:
                        continue
                    close = z1[row] - 1 - frame0
                    fitted = max(close, deadlines[cell] - 1 - frame0)
                    usage[close] += transform
                    usage[fitted + 1] -= transform
            current = peak = 0
            for u in usage:
                current += u
                peak = max(peak, current)
            required = max(required, peak)
        return required

    def task(self, index):
        """
        Get a task.
        """
        job = self.one_pass_jobs[index]
        reflections = self.reflections.select(job.rows)
        reflections[self.LEARN_ONLY] = job.learn_only
        return Task(
            index=index,
            job=job.frames,
            experiments=self.experiments,
            reflections=reflections,
            params=self.params,
            executor=self.executor,
        )

    def accumulate(self, result):
        """
        Accumulate the results: the rows each job owns.
        """
        job = self.one_pass_jobs[result.index]
        assert not self.finished_jobs[result.index]
        reflections = result.reflections
        owned = ~reflections[self.LEARN_ONLY]
        del reflections[self.LEARN_ONLY]
        self.reflections.set_selected(job.rows.select(owned), reflections.select(owned))
        self.data[result.index] = result.data
        self.finished_jobs[result.index] = True
        self.time.read += result.read_time
        self.time.extract += result.extract_time
        self.time.process += result.process_time
        self.time.total += result.total_time

    def finalize(self):
        """
        Finalize the processing and finish.
        """
        start_time = time()
        assert all(self.finished_jobs), "Manager is not finished"
        self.time.finalize = time() - start_time
        self.finalized = True

    def result(self):
        """
        Return the result.
        """
        assert self.finalized, "Manager is not finalized"
        return self.reflections, self.data

    def finished(self):
        return self.finalized and all(self.finished_jobs)

    def __len__(self):
        return len(self.one_pass_jobs)

    def summary(self):
        """
        Get a summary of the processing
        """
        rows = [
            [
                "#",
                "Group",
                "Frame From",
                "Frame To",
                "Angle From",
                "Angle To",
                "# Reflections",
                "# Reference only",
            ]
        ]
        for i, job in enumerate(self.one_pass_jobs):
            f0, f1 = job.frames
            scan = self.experiments[job.expr[0]].scan
            p0 = scan.get_angle_from_array_index(f0)
            p1 = scan.get_angle_from_array_index(f1)
            n_learn = job.learn_only.count(True)
            rows.append(
                [
                    str(i),
                    str(job.group),
                    str(f0 + 1),
                    str(f1),
                    str(p0),
                    str(p1),
                    str(len(job.rows) - n_learn),
                    str(n_learn),
                ]
            )
        return (
            "Processing reflections in one pass, in the following blocks of images:\n\n"
            "{}\n"
        ).format(tabulate(rows, headers="firstrow"))


class OnePassProcessor3D(_ProcessorRot):
    """Top level processor for 3D processing in one pass over the images."""

    def __init__(self, experiments, reflections, params, profile_modeller):
        """Initialise the manager and the processor."""
        params.shoebox.partials = False
        params.shoebox.flatten = False
        manager = _OnePassManager(experiments, reflections, params, profile_modeller)
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
