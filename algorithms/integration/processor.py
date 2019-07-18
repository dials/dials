#
# processor.py
#
#  Copyright (C) 2015 Diamond Light Source
#
#  Author: James Parkhurst
#
#  This code is distributed under the BSD license, a copy of which is
#  included in the root directory of this package.

from __future__ import absolute_import, division, print_function

import logging
import math
from time import time

import boost.python
import libtbx
from dials_algorithms_integration_integrator_ext import *

logger = logging.getLogger(__name__)


class ExecutorAux(Executor, boost.python.injector):
    def __getinitargs__(self):
        return ()


class _Job(object):
    def __init__(self):
        self.index = 0
        self.nthreads = 1


job = _Job()


class MultiProcessing(object):
    """
    Multi processing parameters

    """

    def __init__(self):
        self.method = "multiprocessing"
        self.nproc = 1
        self.njobs = 1
        self.nthreads = 1

    def update(self, other):
        self.method = other.method
        self.nproc = other.nproc
        self.njobs = other.njobs
        self.nthreads = other.nthreads


class Lookup(object):
    """
    Lookup parameters

    """

    def __init__(self):
        self.mask = None

    def update(self, other):
        self.mask = other.mask


class Block(object):
    """
    Block parameters

    """

    def __init__(self):
        self.size = libtbx.Auto
        self.units = "degrees"
        self.threshold = 0.99
        self.force = False
        self.max_memory_usage = 0.75

    def update(self, other):
        self.size = other.size
        self.units = other.units
        self.threshold = other.threshold
        self.force = other.force
        self.max_memory_usage = other.max_memory_usage


class Shoebox(object):
    """
    Shoebox parameters

    """

    def __init__(self):
        self.flatten = False
        self.partials = False

    def update(self, other):
        self.flatten = other.flatten
        self.partials = other.partials


class Debug(object):
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


class Parameters(object):
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


class TimingInfo(object):
    """
    A class to contain timing info.

    """

    def __init__(self):
        self.read = 0
        self.extract = 0
        self.initialize = 0
        self.process = 0
        self.finalize = 0
        self.total = 0
        self.user = 0

    def __str__(self):
        """ Convert to string. """
        from libtbx.table_utils import format as table

        rows = [
            ["Read time", "%.2f seconds" % (self.read)],
            ["Extract time", "%.2f seconds" % (self.extract)],
            ["Pre-process time", "%.2f seconds" % (self.initialize)],
            ["Process time", "%.2f seconds" % (self.process)],
            ["Post-process time", "%.2f seconds" % (self.finalize)],
            ["Total time", "%.2f seconds" % (self.total)],
            ["User time", "%.2f seconds" % (self.user)],
        ]
        return table(rows, justify="right", prefix=" ")


class ExecuteParallelTask(object):
    """
    Helper class to run things on cluster

    """

    def __call__(self, task):
        from dials.util import log

        log.config_simple_cached()
        result = task()
        handlers = logging.getLogger("dials").handlers
        assert len(handlers) == 1, "Invalid number of logging handlers"
        return result, handlers[0].messages()


class Processor(object):
    """ Processor interface class. """

    def __init__(self, manager):
        """
        Initialise the processor.

        The processor requires a manager class implementing the Manager interface.
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
        from dials.util.mp import multi_node_parallel_map
        import platform

        start_time = time()
        self.manager.initialize()
        mp_method = self.manager.params.mp.method
        mp_njobs = self.manager.params.mp.njobs
        mp_nproc = self.manager.params.mp.nproc
        if (
            mp_njobs * mp_nproc
        ) > 1 and platform.system() == "Windows":  # platform.system() forks which is bad for MPI, so don't use it unless nproc > 1
            logger.warning(
                "\n"
                + "*" * 80
                + "\n"
                + "Multiprocessing is not available on windows. Setting nproc = 1\n"
                + "*" * 80
                + "\n"
            )
            mp_nproc = 1
            mp_njobs = 1
        assert mp_nproc > 0, "Invalid number of processors"
        if mp_nproc * mp_njobs > len(self.manager):
            mp_nproc = min(mp_nproc, len(self.manager))
            mp_njobs = int(math.ceil(len(self.manager) / mp_nproc))
        logger.info(self.manager.summary())
        if mp_njobs > 1:
            assert mp_method != "none" and mp_method is not None
            logger.info(
                " Using %s with %d parallel job(s) and %d processes per node\n"
                % (mp_method, mp_njobs, mp_nproc)
            )
        else:
            logger.info(" Using multiprocessing with %d parallel job(s)\n" % (mp_nproc))
        if mp_njobs * mp_nproc > 1:

            def process_output(result):
                for message in result[1]:
                    logger.log(message.levelno, message.msg)
                self.manager.accumulate(result[0])
                result[0].reflections = None
                result[0].data = None

            multi_node_parallel_map(
                func=ExecuteParallelTask(),
                iterable=list(self.manager.tasks()),
                njobs=mp_njobs,
                nproc=mp_nproc,
                callback=process_output,
                cluster_method=mp_method,
                preserve_order=True,
                preserve_exception_message=True,
            )
        else:
            for task in self.manager.tasks():
                self.manager.accumulate(task())
        self.manager.finalize()
        end_time = time()
        self.manager.time.user_time = end_time - start_time
        result1, result2 = self.manager.result()
        return result1, result2, self.manager.time


class Result(object):
    """
    A class representing a processing result.

    """

    def __init__(self, index, reflections, data=None):
        """
        Initialise the data.

        :param index: The processing job index
        :param reflections: The processed reflections
        :param data: Other processed data

        """
        self.index = index
        self.reflections = reflections
        self.data = data


class NullTask(object):
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
        result = Result(self.index, self.reflections, None)
        result.read_time = 0
        result.extract_time = 0
        result.process_time = 0
        result.total_time = 0
        return result


class Task(object):
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
        :param save_shoeboxes: Save the shoeboxes to file
        :param executor: The executor class

        """
        assert executor is not None, "No executor given"
        assert len(reflections) > 0, "Zero reflections given"
        assert params.block.max_memory_usage > 0.0, "Max memory % must be > 0"
        assert params.block.max_memory_usage <= 1.0, "Max memory % must be < 1"
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
        from dials.array_family import flex
        from dials.model.data import make_image
        from libtbx.introspection import machine_memory_info

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
        frame00, frame01 = self.job
        try:
            frame10, frame11 = imageset.get_array_range()
        except Exception:
            frame10, frame11 = (0, len(imageset))
        try:
            assert frame00 < frame01
            assert frame10 < frame11
            assert frame00 >= frame10
            assert frame01 <= frame11
            index0 = frame00 - frame10
            index1 = index0 + (frame01 - frame00)
            assert index0 < index1
            assert index0 >= 0
            assert index1 <= len(imageset)
            imageset = imageset[index0:index1]
        except Exception:
            raise RuntimeError("Programmer Error: bad array range")
        try:
            frame0, frame1 = imageset.get_array_range()
        except Exception:
            frame0, frame1 = (0, len(imageset))

        # Initlize the executor
        self.executor.initialize(frame0, frame1, self.reflections)

        # Set the shoeboxes (dont't allocate)
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

        # Compute percentage of max available. The function is not portable to
        # windows so need to add a check if the function fails. On windows no
        # warning will be printed
        memory_info = machine_memory_info()
        total_memory = memory_info.memory_total()
        sbox_memory = processor.compute_max_memory_usage()
        if total_memory is not None:
            assert total_memory > 0, "Your system appears to have no memory!"
            assert (
                self.params.block.max_memory_usage > 0.0
            ), "maximum memory usage must be > 0"
            assert (
                self.params.block.max_memory_usage <= 1.0
            ), "maximum memory usage must be <= 1"
            limit_memory = total_memory * self.params.block.max_memory_usage
            if sbox_memory > limit_memory:
                raise RuntimeError(
                    """
        There was a problem allocating memory for shoeboxes. Possible solutions
        include increasing the percentage of memory allowed for shoeboxes or
        decreasing the block size. This could also be caused by a highly mosaic
        crystal model - is your crystal really this mosaic?
          Total system memory: %g GB
          Limit shoebox memory: %g GB
          Required shoebox memory: %g GB
        """
                    % (total_memory / 1e9, limit_memory / 1e9, sbox_memory / 1e9)
                )
            else:
                logger.info(" Memory usage:")
                logger.info("  Total system memory: %g GB" % (total_memory / 1e9))
                logger.info("  Limit shoebox memory: %g GB" % (limit_memory / 1e9))
                logger.info("  Required shoebox memory: %g GB" % (sbox_memory / 1e9))
                logger.info("")

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
                        % (len(mask), len(self.params.lookup.mask))
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
                    table.as_pickle("shoeboxes_%d_%d.refl" % (self.index, i))
            else:
                output.as_pickle("shoeboxes_%d.refl" % self.index)

        # Delete the shoeboxes
        if self.params.debug.separate_files or not self.params.debug.output:
            del self.reflections["shoebox"]

        # Finalize the executor
        self.executor.finalize()

        # Return the result
        result = Result(self.index, self.reflections, self.executor.data())
        result.read_time = read_time
        result.extract_time = processor.extract_time()
        result.process_time = processor.process_time()
        result.total_time = time() - start_time
        return result


class Manager(object):
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
        self.time = TimingInfo()

    def initialize(self):
        """
        Initialise the processing

        """
        # Get the start time
        start_time = time()

        # Ensure the reflections contain bounding boxes
        assert "bbox" in self.reflections, "Reflections have no bbox"

        # Compute the block size and processors
        self.compute_blocks()
        self.compute_jobs()
        self.split_reflections()
        self.compute_processors()

        # Create the reflection manager
        self.manager = ReflectionManager(self.jobs, self.reflections)

        # Parallel reading of HDF5 from the same handle is not allowed. Python
        # multiprocessing is a bit messed up and used fork on linux so need to
        # close and reopen file.
        for exp in self.experiments:
            if exp.imageset.reader().is_single_file_reader():
                exp.imageset.reader().nullify_format_instance()

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
            logger.warning("*** WARNING: no reflections in job %d ***", index)
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
        """ Accumulate the results. """
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

    def compute_blocks(self):
        """
        Compute the processing block size.

        """

        if self.params.block.size == libtbx.Auto:
            if (
                self.params.mp.nproc * self.params.mp.njobs == 1
                and not self.params.debug.output
                and not self.params.block.force
            ):
                self.params.block.size = None
            else:
                assert self.params.block.threshold > 0, "Threshold must be > 0"
                assert self.params.block.threshold <= 1.0, "Threshold must be < 1"
                nframes = sorted([b[5] - b[4] for b in self.reflections["bbox"]])
                cutoff = int(self.params.block.threshold * len(nframes))
                block_size = nframes[cutoff] * 2
                self.params.block.size = block_size
                self.params.block.units = "frames"

    def compute_jobs(self):
        """
        Compute the jobs

        """
        from itertools import groupby

        groups = groupby(
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
            elif self.params.block.units == "radians":
                phi0, dphi = scan.get_oscillation(deg=False)
                block_size_frames = int(math.ceil(self.params.block.size / dphi))
            elif self.params.block.units == "degrees":
                phi0, dphi = scan.get_oscillation()
                block_size_frames = int(math.ceil(self.params.block.size / dphi))
            elif self.params.block.units == "frames":
                block_size_frames = int(math.ceil(self.params.block.size))
            else:
                raise RuntimeError("Unknown block_size_units = %s" % block_size_units)
            self.jobs.add((i0, i1), array_range, block_size_frames)
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
                    " Split %d reflections into %d partial reflections\n"
                    % (num_full, num_partial)
                )
        else:
            num_full = len(self.reflections)
            self.jobs.split(self.reflections)
            num_partial = len(self.reflections)
            assert num_partial >= num_full, "Invalid number of partials"
            if num_partial > num_full:
                num_split = num_partial - num_full
                logger.info(
                    " Split %d reflections overlapping job boundaries\n" % num_split
                )

        # Compute the partiality
        self.reflections.compute_partiality(self.experiments)

    def compute_processors(self):
        """
        Compute the number of processors

        """
        from libtbx.introspection import machine_memory_info
        from dials.array_family import flex

        # Set the memory usage per processor
        if self.params.mp.method == "multiprocessing" and self.params.mp.nproc > 1:

            # Get the maximum shoebox memory
            max_memory = flex.max(
                self.jobs.shoebox_memory(self.reflections, self.params.shoebox.flatten)
            )

            # Compute percentage of max available. The function is not portable to
            # windows so need to add a check if the function fails. On windows no
            # warning will be printed
            memory_info = machine_memory_info()
            total_memory = memory_info.memory_total()
            if total_memory is not None:
                assert total_memory > 0, "Your system appears to have no memory!"
                limit_memory = total_memory * self.params.block.max_memory_usage
                njobs = int(math.floor(limit_memory / max_memory))
                if njobs < 1:
                    raise RuntimeError(
                        """
            No enough memory to run integration jobs. Possible solutions
            include increasing the percentage of memory allowed for shoeboxes or
            decreasing the block size.
              Total system memory: %g GB
              Limit shoebox memory: %g GB
              Max shoebox memory: %g GB
          """
                        % (total_memory / 1e9, limit_memory / 1e9, max_memory / 1e9)
                    )
                else:
                    self.params.mp.nproc = min(self.params.mp.nproc, njobs)
                    self.params.block.max_memory_usage /= self.params.mp.nproc

    def summary(self):
        """
        Get a summary of the processing

        """
        from libtbx.table_utils import format as table

        # Compute the task table
        if self.experiments.all_stills():
            rows = [["#", "Group", "Frame From", "Frame To", "# Reflections"]]
            for i in range(len(self)):
                job = self.manager.job(i)
                group = job.index()
                f0, f1 = job.frames()
                n = self.manager.num_reflections(i)
                rows.append([str(i), str(group), str(f0), str(f1), str(n)])
        elif self.experiments.all_sweeps():
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
                    [str(i), str(group), str(f0), str(f1), str(p0), str(p1), str(n)]
                )
        else:
            raise RuntimeError("Experiments must be all sweeps or all stills")

        # The job table
        task_table = table(rows, has_header=True, justify="right", prefix=" ")

        # The format string
        if self.params.block.size is None:
            block_size = "auto"
        else:
            block_size = str(self.params.block.size)
        fmt = (
            "Processing reflections in the following blocks of images:\n"
            "\n"
            " block_size: %s %s\n"
            "\n"
            "%s\n"
        )
        return fmt % (block_size, self.params.block.units, task_table)


class ManagerRot(Manager):
    """Specialize the manager for oscillation data using the oscillation pre and
    post processors."""

    def __init__(self, experiments, reflections, params):
        """ Initialise the pre-processor, post-processor and manager. """

        # Ensure we have the correct type of data
        if not experiments.all_sweeps():
            raise RuntimeError(
                """
        An inappropriate processing algorithm may have been selected!
         Trying to perform rotation processing when not all experiments
         are indicated as rotation experiments.
      """
            )

        # Initialise the manager
        super(ManagerRot, self).__init__(experiments, reflections, params)


class ManagerStills(Manager):
    """Specialize the manager for stills data using the stills pre and post
    processors."""

    def __init__(self, experiments, reflections, params):
        """ Initialise the pre-processor, post-processor and manager. """

        # Ensure we have the correct type of data
        if not experiments.all_stills():
            raise RuntimeError(
                """
        An inappropriate processing algorithm may have been selected!
         Trying to perform stills processing when not all experiments
         are indicated as stills experiments.
      """
            )

        # Initialise the manager
        super(ManagerStills, self).__init__(experiments, reflections, params)


class Processor3D(Processor):
    """ Top level processor for 3D processing. """

    def __init__(self, experiments, reflections, params):
        """ Initialise the manager and the processor. """

        # Set some parameters
        params.shoebox.partials = False
        params.shoebox.flatten = False

        # Create the processing manager
        manager = ManagerRot(experiments, reflections, params)

        # Initialise the processor
        super(Processor3D, self).__init__(manager)


class ProcessorFlat3D(Processor):
    """ Top level processor for flat 2D processing. """

    def __init__(self, experiments, reflections, params):
        """ Initialise the manager and the processor. """

        # Set some parameters
        params.shoebox.flatten = True
        params.shoebox.partials = False

        # Create the processing manager
        manager = ManagerRot(experiments, reflections, params)

        # Initialise the processor
        super(ProcessorFlat3D, self).__init__(manager)


class Processor2D(Processor):
    """ Top level processor for 2D processing. """

    def __init__(self, experiments, reflections, params):
        """ Initialise the manager and the processor. """

        # Set some parameters
        params.shoebox.partials = True

        # Create the processing manager
        manager = ManagerRot(experiments, reflections, params)

        # Initialise the processor
        super(Processor2D, self).__init__(manager)


class ProcessorSingle2D(Processor):
    """ Top level processor for still image processing. """

    def __init__(self, experiments, reflections, params):
        """ Initialise the manager and the processor. """

        # Set some of the parameters
        params.block.size = 1
        params.block.units = "frames"
        params.shoebox.partials = True
        params.shoebox.flatten = False

        # Create the processing manager
        manager = ManagerRot(experiments, reflections, params)

        # Initialise the processor
        super(ProcessorSingle2D, self).__init__(manager)


class ProcessorStills(Processor):
    """ Top level processor for still image processing. """

    def __init__(self, experiments, reflections, params):
        """ Initialise the manager and the processor. """

        # Set some parameters
        params.block.size = 1
        params.block.units = "frames"
        params.shoebox.partials = False
        params.shoebox.flatten = False

        # Create the processing manager
        manager = ManagerStills(experiments, reflections, params)

        # Initialise the processor
        super(ProcessorStills, self).__init__(manager)


class ProcessorBuilder(object):
    """
    A class to simplify building the processor

    """

    def __init__(self, Class, experiments, reflections, params=None):
        """
        Initialize with the required input

        :param Class: The input class
        :param experiments: The input experiments
        :param reflections: The reflections
        :param params: Optionally input parameters

        """
        self.Class = Class
        self.experiments = experiments
        self.reflections = reflections
        self.params = Parameters()
        if params is not None:
            self.params.update(params)

    def build(self):
        """
        Build the class

        :return: The processor class

        """
        return self.Class(self.experiments, self.reflections, self.params)
