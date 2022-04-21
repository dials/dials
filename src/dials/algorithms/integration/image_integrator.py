from __future__ import annotations

import logging
import platform
from time import time

import dials.algorithms.integration
from dials.algorithms.integration.processor import job
from dials.model.data import ImageVolume, MultiPanelImageVolume, make_image
from dials.util import log
from dials.util.log import rehandle_cached_records
from dials.util.mp import multi_node_parallel_map
from dials_algorithms_integration_integrator_ext import ReflectionManagerPerImage

logger = logging.getLogger(__name__)


class ProcessorImage:
    """Top level processor for per image processing."""

    def __init__(self, experiments, reflections, params):
        """
        Initialise the manager and the processor.

        The processor requires a manager class implementing the Manager interface.
        This class executes all the workers in separate threads and accumulates the
        results to expose to the user.

        :param params: The phil parameters
        """

        # Create the processing manager
        self.manager = ManagerImage(experiments, reflections, params)

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
        mp_method = self.manager.params.integration.mp.method
        mp_nproc = min(len(self.manager), self.manager.params.integration.mp.nproc)
        mp_njobs = self.manager.params.integration.mp.njobs
        if (
            mp_nproc > 1 and platform.system() == "Windows"
        ):  # platform.system() forks which is bad for MPI, so don't use it unless nproc > 1
            logger.warning(
                "Multiprocessing is not available on windows. Setting nproc = 1\n"
            )
            mp_nproc = 1
        assert mp_nproc > 0, "Invalid number of processors"
        logger.info(self.manager.summary())
        logger.info(" Using %s with %d parallel job(s)\n", mp_method, mp_nproc)
        if mp_nproc > 1:

            def process_output(result):
                rehandle_cached_records(result[1])
                self.manager.accumulate(result[0])
                result[0].reflections = None
                result[0].data = None

            def execute_task(task):
                log.config_simple_cached()
                result = task()
                handlers = logging.getLogger("dials").handlers
                assert len(handlers) == 1, "Invalid number of logging handlers"
                return result, handlers[0].records

            multi_node_parallel_map(
                func=execute_task,
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
        result = self.manager.result()
        return result, self.manager.time


class Task:
    """
    A class to perform a null task.
    """

    def __init__(self, index, frames, reflections, experiments, params, executor):
        """
        Initialise the task

        :param index: The index of the processing job
        :param frames: The frames to process
        :param experiments: The list of experiments
        :param reflections: The list of reflections
        :param params The processing parameters
        :param executor: The executor class
        """
        self.index = index
        self.frames = frames
        self.experiments = experiments
        self.reflections = reflections
        self.params = params
        self.executor = executor

    def __call__(self):
        """
        Do the processing.

        :return: The processed data
        """
        # Set the job index
        job.index = self.index

        # Get the start time
        start_time = time()

        # Check all reflections have same imageset and get it
        exp_id = list(set(self.reflections["id"]))
        imageset = self.experiments[exp_id[0]].imageset
        for i in exp_id[1:]:
            assert (
                self.experiments[i].imageset == imageset
            ), "Task can only handle 1 imageset"

        # Get the sub imageset
        frame00, frame01 = self.frames
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

        # Initialise the dataset
        image_volume = MultiPanelImageVolume()
        for panel in self.experiments[0].detector:
            image_volume.add(
                ImageVolume(
                    frame0, frame1, panel.get_image_size()[1], panel.get_image_size()[0]
                )
            )

        # Read all the images into a block of data
        read_time = 0.0
        for i in range(len(imageset)):
            st = time()
            image = imageset.get_corrected_data(i)
            mask = imageset.get_mask(i)
            if self.params.integration.lookup.mask is not None:
                assert len(mask) == len(
                    self.params.lookup.mask
                ), "Mask/Image are incorrect size %d %d" % (
                    len(mask),
                    len(self.params.integration.lookup.mask),
                )
                mask = tuple(
                    m1 & m2 for m1, m2 in zip(self.params.integration.lookup.mask, mask)
                )
            image_volume.set_image(frame0 + i, make_image(image, mask))
            read_time += time() - st
            del image
            del mask

        # Process the data
        st = time()
        data = self.executor.process(image_volume, self.experiments, self.reflections)
        process_time = time() - st

        # Set the result values
        return dials.algorithms.integration.Result(
            index=self.index,
            reflections=self.reflections,
            read_time=read_time,
            process_time=process_time,
            total_time=time() - start_time,
            extract_time=0,
            data=data,
        )


class ManagerImage:
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

        # Split the reflections into partials
        self._split_reflections()

        # Create the reflection manager
        frames = self.experiments[0].scan.get_array_range()
        self.manager = ReflectionManagerPerImage(frames, self.reflections)

        # Set the initialization time
        self.time.initialize = time() - start_time

    def task(self, index):
        """
        Get a task.
        """
        return Task(
            index=index,
            frames=self.manager.frames(index),
            reflections=self.manager.split(index),
            experiments=self.experiments,
            params=self.params,
            executor=self.executor,
        )

    def tasks(self):
        """
        Iterate through the tasks.
        """
        for i in range(len(self)):
            yield self.task(i)

    def accumulate(self, result):
        """
        Accumulate the results.
        """
        self.manager.accumulate(result.index, result.reflections)
        if result.data is not None:
            self.executor.accumulate(result.index, result.data)
        self.time.read += result.read_time
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
        return self.reflections

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

    def summary(self):
        return ""

    def _split_reflections(self):
        """
        Split the reflections into partials or over job boundaries
        """

        # Optionally split the reflection table into partials, otherwise,
        # split over job boundaries
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
