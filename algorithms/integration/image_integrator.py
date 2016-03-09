#
# image_integrator.py
#
#  Copyright (C) 2013 Diamond Light Source
#
#  Author: James Parkhurst
#
#  This code is distributed under the BSD license, a copy of which is
#  included in the root directory of this package.


from __future__ import division


class TimingInfo(object):
  '''
  A class to contain timing info.

  '''
  def __init__(self):
    self.read = 0
    self.initialize = 0
    self.process = 0
    self.finalize = 0
    self.total = 0
    self.user = 0

  def __str__(self):
    ''' Convert to string. '''
    from libtbx.table_utils import format as table
    rows = [
      ["Read time"        , "%.2f seconds" % (self.read)       ],
      ["Pre-process time" , "%.2f seconds" % (self.initialize) ],
      ["Process time"     , "%.2f seconds" % (self.process)    ],
      ["Post-process time", "%.2f seconds" % (self.finalize)   ],
      ["Total time"       , "%.2f seconds" % (self.total)      ],
      ["User time"        , "%.2f seconds" % (self.user)       ],
    ]
    return table(rows, justify='right', prefix=' ')


class ProcessorImageBase(object):
  ''' Processor interface class. '''

  def __init__(self, manager):
    '''
    Initialise the processor.

    The processor requires a manager class implementing the Manager interface.
    This class executes all the workers in separate threads and accumulates the
    results to expose to the user.

    :param manager: The processing manager
    :param params: The phil parameters

    '''
    self.manager = manager

  @property
  def executor(self):
    '''
    Get the executor

    :return: The executor

    '''
    return self.manager.executor

  @executor.setter
  def executor(self, function):
    '''
    Set the executor

    :param function: The executor

    '''
    self.manager.executor = function

  def process(self):
    '''
    Do all the processing tasks.

    :return: The processing results

    '''
    from time import time
    from libtbx import easy_mp
    from logging import info, warn
    import platform
    start_time = time()
    self.manager.initialize()
    mp_method = self.manager.params.integration.mp.method
    mp_nproc = min(len(self.manager), self.manager.params.integration.mp.nproc)
    if platform.system() == "Windows" and mp_nproc > 1:
      warn("")
      warn("*" * 80)
      warn("Multiprocessing is not available on windows. Setting nproc = 1")
      warn("*" * 80)
      warn("")
      mp_nproc = 1
    assert mp_nproc > 0, "Invalid number of processors"
    info(self.manager.summary())
    info(' Using %s with %d parallel job(s)\n' % (
      mp_method, mp_nproc))
    if mp_nproc > 1:
      def process_output(result):
        import logging
        for message in result[1]:
          logging.log(message.levelno, message.msg)
        self.manager.accumulate(result[0])
        result[0].reflections = None
        result[0].data = None
      def execute_task(task):
        from cStringIO import StringIO
        from dials.util import log
        import logging
        log.config_simple_cached()
        result = task()
        handlers = logging.getLogger().handlers
        assert len(handlers) == 1, "Invalid number of logging handlers"
        return result, handlers[0].messages()
      easy_mp.parallel_map(
        func=execute_task,
        iterable=list(self.manager.tasks()),
        processes=mp_nproc,
        callback=process_output,
        method=mp_method,
        preserve_order=True,
        preserve_exception_message=True)
    else:
      for task in self.manager.tasks():
        self.manager.accumulate(task())
    self.manager.finalize()
    end_time = time()
    self.manager.time.user_time = end_time - start_time
    result = self.manager.result()
    return result, self.manager.time


class Result(object):
  '''
  A class representing a processing result.

  '''

  def __init__(self, index, reflections):
    '''
    Initialise the data.

    :param index: The processing job index
    :param reflections: The processed reflections
    :param data: Other processed data

    '''
    self.index = index
    self.reflections = reflections


class Dataset(object):

  def __init__(self, nframes, size):
    from dials.array_family import flex
    self.data = []
    self.mask = []
    for sz in size:
      self.data.append(flex.double(flex.grid(nframes, sz[0], sz[1])))
      self.mask.append(flex.bool(flex.grid(nframes, sz[0], sz[1])))

  def set_image(self, index, data, mask):
    from dials.array_family import flex
    for d1, d2 in zip(self.data, data):
      h,w = d2.all()
      d2.reshape(flex.grid(1,h,w))
      d1[index:index+1,:,:] = d2.as_double()
    for m1, m2 in zip(self.mask, mask):
      h,w = m2.all()
      m2.reshape(flex.grid(1,h,w))
      m1[index:index+1,:,:] = m2

class Task(object):
  '''
  A class to perform a null task.

  '''
  def __init__(self,
               index,
               frames,
               reflections,
               experiments,
               params,
               executor):
    '''
    Initialise the task

    :param index: The index of the processing job
    :param frames: The frames to process
    :param experiments: The list of experiments
    :param reflections: The list of reflections
    :param params The processing parameters
    :param executor: The executor class

    '''
    self.index = index
    self.frames = frames
    self.experiments = experiments
    self.reflections = reflections
    self.params = params
    self.executor = executor

  def __call__(self):
    '''
    Do the processing.

    :return: The processed data

    '''
    from time import time
    from logging import info

    # Get the start time
    start_time = time()

    # Check all reflections have same imageset and get it
    exp_id = list(set(self.reflections['id']))
    imageset = self.experiments[exp_id[0]].imageset
    for i in exp_id[1:]:
      assert self.experiments[i].imageset == imageset, "Task can only handle 1 imageset"

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
      raise RuntimeError('Programmer Error: bad array range')
    try:
      frame0, frame1 = imageset.get_array_range()
    except Exception:
      frame0, frame1 = (0, len(imageset))

    # Initialise the executor
    self.executor.initialize()

    # Initialise the dataset
    size = []
    for panel in self.experiments[0].detector:
      size.append(panel.get_image_size()[::-1])
    dataset = Dataset(len(imageset), size)

    # Read all the images into a block of data
    read_time = 0.0
    for i in range(len(imageset)):
      st = time()
      info("Reading image %d" % (self.frames[0] + i))
      image = imageset.get_corrected_data(i)
      mask  = imageset.get_mask(i)
      if self.params.integration.lookup.mask is not None:
        assert len(mask) == len(self.params.lookup.mask), \
          "Mask/Image are incorrect size %d %d" % (
            len(mask),
            len(self.params.lookup.mask))
        mask = tuple(m1 & m2 for m1, m2 in zip(self.params.lookup.mask, mask))
      dataset.set_image(i, image, mask)
      read_time += time() - st
      del image
      del mask

    # Process the data
    st = time()
    self.executor.process(dataset)
    process_time = time() - st

    # Finalize the executor
    self.executor.finalize()

    # Set the result values
    result = Result(self.index, self.reflections)
    result.read_time = read_time
    result.process_time = process_time
    result.total_time = time() - start_time
    return result


class ManagerImage(object):
  '''
  A class to manage processing book-keeping

  '''

  def __init__(self, experiments, reflections, params):
    '''
    Initialise the manager.

    :param experiments: The list of experiments
    :param reflections: The list of reflections
    :param params: The phil parameters

    '''
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
    self.time = TimingInfo()

  def initialize(self):
    '''
    Initialise the processing

    '''
    from dials_algorithms_integration_integrator_ext import ReflectionManagerPerImage
    from time import time

    # Get the start time
    start_time = time()

    # Ensure the reflections contain bounding boxes
    assert "bbox" in self.reflections, "Reflections have no bbox"

    # Create the reflection manager
    frames = self.experiments[0].scan.get_array_range()
    self.manager = ReflectionManagerPerImage(frames, self.reflections)

    # Parallel reading of HDF5 from the same handle is not allowed. Python
    # multiprocessing is a bit messed up and used fork on linux so need to
    # close and reopen file.
    from dxtbx.imageset import SingleFileReader
    for exp in self.experiments:
      if isinstance(exp.imageset.reader(), SingleFileReader):
        exp.imageset.reader().nullify_format_instance()

    # Set the initialization time
    self.time.initialize = time() - start_time

  def task(self, index):
    '''
    Get a task.

    '''
    return Task(
        index       = index,
        frames      = self.manager.frames(index),
        reflections = self.manager.split(index),
        experiments = self.experiments,
        params      = self.params,
        executor    = self.executor)

  def tasks(self):
    '''
    Iterate through the tasks.

    '''
    for i in range(len(self)):
      yield self.task(i)

  def accumulate(self, result):
    '''
    Accumulate the results.

    '''
    self.manager.accumulate(result.index, result.reflections)
    self.time.read += result.read_time
    self.time.process += result.process_time
    self.time.total += result.total_time

  def finalize(self):
    '''
    Finalize the processing and finish.

    '''
    from time import time

    # Get the start time
    start_time = time()

    # Check manager is finished
    assert self.manager.finished(), "Manager is not finished"

    # Update the time and finalized flag
    self.time.finalize = time() - start_time
    self.finalized = True

  def result(self):
    '''
    Return the result.

    :return: The result

    '''
    assert self.finalized == True, "Manager is not finalized"
    return self.reflections

  def finished(self):
    '''
    Return if all tasks have finished.

    :return: True/False all tasks have finished

    '''
    return self.finalized and self.manager.finished()

  def __len__(self):
    '''
    Return the number of tasks.

    :return: the number of tasks

    '''
    return len(self.manager)

  def summary(self):
    return ''


class ProcessorImage(ProcessorImageBase):
  ''' Top level processor for per image processing. '''

  def __init__(self, experiments, reflections, params):
    ''' Initialise the manager and the processor. '''

    # Create the processing manager
    manager = ManagerImage(experiments, reflections, params)

    # Initialise the processor
    super(ProcessorImage, self).__init__(manager)


class ImageIntegratorExecutor(object):

  def __init__(self):
    pass

  def initialize(self):
    pass

  def process(self, dataset):
    pass

  def finalize(self):
    pass


class ImageIntegrator(object):
  '''
  A class that does integration directly on the image skipping the shoebox
  creation step.

  '''

  def __init__(self,
               experiments,
               reflections,
               params):
    '''
    Initialize the integrator

    :param experiments: The experiment list
    :param reflections: The reflections to process
    :param params: The parameters to use

    '''

    # Save some stuff
    self.experiments = experiments
    self.reflections = reflections
    self.params = params
    self.profile_model_report = None
    self.integration_report = None

  def integrate(self):
    '''
    Integrate the data

    '''
    from dials.algorithms.integration.report import IntegrationReport
    from dials.util.command_line import heading
    from logging import info, debug

    # Init the report
    self.profile_model_report = None
    self.integration_report = None

    # Heading
    info("=" * 80)
    info("")
    info(heading("Processing reflections"))
    info("")

    # Create summary format
    fmt = (
      ' Processing the following experiments:\n'
      '\n'
      ' Experiments: %d\n'
      ' Beams:       %d\n'
      ' Detectors:   %d\n'
      ' Goniometers: %d\n'
      ' Scans:       %d\n'
      ' Crystals:    %d\n'
      ' Imagesets:   %d\n'
    )

    # Print the summary
    info(fmt % (
      len(self.experiments),
      len(self.experiments.beams()),
      len(self.experiments.detectors()),
      len(self.experiments.goniometers()),
      len(self.experiments.scans()),
      len(self.experiments.crystals()),
      len(self.experiments.imagesets())))

    # Print a heading
    info("=" * 80)
    info("")
    info(heading("Integrating reflections"))
    info("")

    # Construvt the image integrator processor
    processor = ProcessorImage(
      self.experiments,
      self.reflections,
      self.params)
    processor.executor = ImageIntegratorExecutor()

    # Do the processing
    self.reflections, time_info = processor.process()

    # Create the integration report
    self.integration_report = IntegrationReport(
      self.experiments,
      self.reflections)
    info("")
    info(self.integration_report.as_str(prefix=' '))

    # Print the time info
    info(str(time_info))
    info("")
