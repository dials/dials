#
# processor.py
#
#  Copyright (C) 2015 Diamond Light Source
#
#  Author: James Parkhurst
#
#  This code is distributed under the BSD license, a copy of which is
#  included in the root directory of this package.

from __future__ import division
from dials_algorithms_integration_integrator_ext import *
from iotbx import phil
import libtbx


class _Job(object):
  def __init__(self):
    self.index = 0
    self.nthreads = 1


job = _Job()


def generate_phil_scope():
  '''
  Generate the processing phil scope.

  :return: The phil scope

  '''

  phil_scope = phil.parse('''

    include scope dials.data.multiprocessing.phil_scope
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

    shoebox {

      flatten = False
        .type = bool
        .help = "Flatten shoeboxes"

      partials = False
        .type = bool
        .help = "Split reflections into partials"

    }

    debug {

      save_shoeboxes = False
        .type = bool
        .help = "Save shoeboxes after each processing task."

      select {

        using = *sum prf
          .type = choice
          .help = "Select using which intensities"

        i_over_sigma_lt = None
          .type = float
          .help = "Select reflections with I/Sigma < value"

        i_over_sigma_gt = None
          .type = float
          .help = "Select reflections with I/Sigma > value"

      }

    }
  ''', process_includes=True)

  # Return the phil scope
  return phil_scope

# The processor phil scope
phil_scope = generate_phil_scope()


class TimingInfo(object):
  '''
  A class to contain timing info.

  '''
  def __init__(self):
    self.read = 0
    self.extract = 0
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
      ["Extract time"     , "%.2f seconds" % (self.extract)    ],
      ["Pre-process time" , "%.2f seconds" % (self.initialize) ],
      ["Process time"     , "%.2f seconds" % (self.process)    ],
      ["Post-process time", "%.2f seconds" % (self.finalize)   ],
      ["Total time"       , "%.2f seconds" % (self.total)      ],
      ["User time"        , "%.2f seconds" % (self.user)       ],
    ]
    return table(rows, justify='right', prefix=' ')

class Processor(object):
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
    from logging import info

    start_time = time()
    self.manager.initialize()
    mp_method = self.manager.mp_method
    mp_nproc = min(len(self.manager), self.manager.mp_nproc)
    mp_nthreads = self.manager.mp_nthreads
    assert mp_nproc > 0, "Invalid number of processors"
    job.nthreads = mp_nthreads
    info(self.manager.summary())
    info(' Using %s with %d parallel job(s) and %d thread(s) per job\n' % (
      mp_method, mp_nproc, mp_nthreads))
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
    result1, result2 = self.manager.result()
    return result1, result2, self.manager.time


class Result(object):
  '''
  A class representing a processing result.

  '''

  def __init__(self, index, reflections, data=None):
    '''
    Initialise the data.

    :param index: The processing job index
    :param reflections: The processed reflections
    :param data: Other processed data

    '''
    self.index = index
    self.reflections = reflections
    self.data = data


class NullTask(object):
  '''
  A class to perform a null task.

  '''
  def __init__(self, index, reflections):
    '''
    Initialise the task

    :param index: The index of the processing job
    :param experiments: The list of experiments
    :param profile_model: The profile model
    :param reflections: The list of reflections

    '''
    self.index = index
    self.reflections = reflections

  def __call__(self):
    '''
    Do the processing.

    :return: The processed data

    '''
    result = Result(self.index, self.reflections, None)
    result.read_time = 0
    result.extract_time = 0
    result.process_time = 0
    result.total_time = 0
    return result


class Task(object):
  '''
  A class to perform a processing task.

  '''

  def __init__(self,
               index,
               job,
               experiments,
               profile_model,
               reflections,
               params,
               mask=None,
               flatten=False,
               save_shoeboxes=False,
               max_memory_usage=0.75,
               executor=None):
    '''
    Initialise the task.

    :param index: The index of the processing job
    :param experiments: The list of experiments
    :param profile_model: The profile model
    :param reflections: The list of reflections
    :param params: The processing parameters
    :param job: The frames to integrate
    :param flatten: Flatten the shoeboxes
    :param save_shoeboxes: Save the shoeboxes to file
    :param executor: The executor class

    '''
    assert executor is not None, "No executor given"
    assert len(reflections) > 0, "Zero reflections given"
    assert max_memory_usage >  0.0, "Max memory % must be > 0"
    assert max_memory_usage <= 1.0, "Max memory % must be < 1"
    self.index = index
    self.job = job
    self.experiments = experiments
    self.profile_model = profile_model
    self.reflections = reflections
    self.params = params
    self.mask = mask
    self.flatten = flatten
    self.save_shoeboxes = save_shoeboxes
    self.max_memory_usage = max_memory_usage
    self.executor = executor

  def __call__(self):
    '''
    Do the processing.

    :return: The processed data

    '''
    from dials.array_family import flex
    from time import time
    from dials.model.data import make_image
    from libtbx.introspection import machine_memory_info
    from logging import info

    # Get the start time
    start_time = time()

    # Set the global process ID
    job.index = self.index

    # Get the sub imageset
    imagesets = self.experiments.imagesets()
    assert len(imagesets) == 1, "Task can only handle 1 imageset"
    imageset = imagesets[0]
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
      raise RuntimeError('Programmer Error: bad array range')
    try:
      frame0, frame1 = imageset.get_array_range()
    except Exception:
      frame0, frame1 = (0, len(imageset))

    # Initlize the executor
    self.executor.initialize(frame0, frame1, self.reflections)

    # Set the shoeboxes (dont't allocate)
    self.reflections['shoebox'] = flex.shoebox(
      self.reflections['panel'],
      self.reflections['bbox'],
      allocate=False,
      flatten=self.flatten)

    # Create the processor
    processor = ShoeboxProcessor(
      self.reflections,
      len(imageset.get_detector()),
      frame0,
      frame1,
      self.save_shoeboxes)

    # Compute percentage of max available. The function is not portable to
    # windows so need to add a check if the function fails. On windows no
    # warning will be printed
    memory_info = machine_memory_info()
    total_memory = memory_info.memory_total()
    sbox_memory = processor.compute_max_memory_usage()
    if total_memory is not None:
      assert total_memory > 0, "Your system appears to have no memory!"
      assert self.max_memory_usage >  0.0
      assert self.max_memory_usage <= 1.0
      limit_memory = total_memory * self.max_memory_usage
      if sbox_memory > limit_memory:
        raise RuntimeError('''
        There was a problem allocating memory for shoeboxes. Possible solutions
        include increasing the percentage of memory allowed for shoeboxes or
        decreasing the block size. This could also be caused by a highly mosaic
        crystal model - is your crystal really this mosaic?
          Total system memory: %g GB
          Limit shoebox memory: %g GB
          Required shoebox memory: %g GB
        ''' % (total_memory/1e9, limit_memory/1e9, sbox_memory/1e9))
      else:
        info(' Memory usage:')
        info('  Total system memory: %g GB' % (total_memory/1e9))
        info('  Limit shoebox memory: %g GB' % (limit_memory/1e9))
        info('  Required shoebox memory: %g GB' % (sbox_memory/1e9))
        info('')

    # Loop through the imageset, extract pixels and process reflections
    read_time = 0.0
    for i in range(len(imageset)):
      st = time()
      image = imageset.get_image(i)
      mask = imageset.get_mask(i)
      if self.mask is not None:
        assert len(mask) == len(self.mask), "Mask/Image are incorrect size"
        mask = tuple(m1 & m2 for m1, m2 in zip(self.mask, mask))
      read_time += time() - st
      processor.next(make_image(image, mask), self.executor)
      del image
      del mask
    assert processor.finished(), "Data processor is not finished"

    # Optionally save the shoeboxes
    if self.save_shoeboxes:
      filename = 'shoeboxes_%d.pickle' % self.index
      output = self.reflections
      if (self.params.debug.select.i_over_sigma_lt is not None or
          self.params.debug.select.i_over_sigma_gt is not None):
        if self.params.debug.select.using == 'sum':
          flag = output.flags.integrated_sum
          Icol = 'intensity.sum.value'
          Vcol = 'intensity.sum.variance'
        else:
          flag = output.flags.integrated_prf
          Icol = 'intensity.prf.value'
          Vcol = 'intensity.prf.variance'
        output = output.select(output.get_flags(flag))
        I = output[Icol]
        V = output[Vcol]
        assert V.all_ge(0), "Some variances < 0"
        output = output.select(V > 0)
        I = output[Icol]
        V = output[Vcol]
        IOS = I / flex.sqrt(V)
        if self.params.debug.select.i_over_sigma_lt is not None:
          mask = IOS < self.params.debug.select.i_over_sigma_lt
          output = output.select(mask)
        if self.params.debug.select.i_over_sigma_gt is not None:
          mask = IOS > self.params.debug.select.i_over_sigma_gt
          output = output.select(mask)
      output.as_pickle(filename)

    # Delete the shoeboxes
    del self.reflections['shoebox']

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
  '''
  A class to manage processing book-keeping

  '''

  def __init__(self, experiments, profile_model, reflections, params):
    '''
    Initialise the manager.

    :param experiments: The list of experiments
    :param profile_model: The profile model
    :param reflections: The list of reflections
    :param params: The phil parameters

    '''

    # Initialise the callbacks
    self.executor = None

    # Save some data
    self.experiments = experiments
    self.profile_model = profile_model
    self.reflections = reflections
    self.mask = params.lookup.mask

    # Other data
    self.data = {}

    # Save some parameters
    self.mp_nproc = params.mp.nproc
    self.partials = params.shoebox.partials
    self.flatten = params.shoebox.flatten
    self.block_size = params.block.size
    self.block_size_units = params.block.units
    self.block_size_threshold = params.block.threshold
    self.block_size_force = params.block.force
    self.save_shoeboxes = params.debug.save_shoeboxes
    self.params = params

    # Save some multiprocessing stuff
    self.mp_nproc = params.mp.nproc
    self.mp_method = params.mp.method
    self.mp_nthreads = params.mp.nthreads

    # Initialise the max memory usage
    self.max_memory_usage = params.block.max_memory_usage

    # Set the finalized flag to False
    self.finalized = False

    # Initialise the timing information
    self.time = TimingInfo()

  def initialize(self):
    '''
    Initialise the processing

    '''
    from time import time

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

    # Set the initialization time
    self.time.initialize = time() - start_time

  def task(self, index):
    '''
    Get a task.

    '''
    from logging import warn
    job = self.manager.job(index)
    frames = job.frames()
    expr_id = job.expr()
    assert expr_id[1] > expr_id[0], "Invalid experiment id"
    assert expr_id[0] >= 0, "Invalid experiment id"
    assert expr_id[1] <= len(self.experiments), "Invalid experiment id"
    expriments = self.experiments[expr_id[0]:expr_id[1]]
    profile_model = self.profile_model[expr_id[0]:expr_id[1]]
    reflections = self.manager.split(index)
    if len(reflections) == 0:
      warn("*** WARNING: no reflections in job %d ***" % index)
      task = NullTask(
        index=index,
        reflections=reflections)
    else:
      task = Task(
        index=index,
        job=frames,
        experiments=expriments,
        profile_model=profile_model,
        reflections=reflections,
        params=self.params,
        mask=self.mask,
        flatten=self.flatten,
        save_shoeboxes=self.save_shoeboxes,
        max_memory_usage=self.max_memory_usage,
        executor=self.executor)
    return task

  def tasks(self):
    '''
    Iterate through the tasks.

    '''
    for i in range(len(self)):
      yield self.task(i)

  def accumulate(self, result):
    ''' Accumulate the results. '''
    self.data[result.index] = result.data
    self.manager.accumulate(result.index, result.reflections)
    self.time.read += result.read_time
    self.time.extract += result.extract_time
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
    return self.manager.data(), self.data

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

  def compute_blocks(self):
    '''
    Compute the processing block size.

    '''
    from logging import info
    from math import ceil, pi
    if self.block_size == libtbx.Auto:
      if (self.mp_nproc == 1 and
          self.save_shoeboxes == False and
          self.block_size_force == False):
        self.block_size = None
      else:
        assert self.block_size_threshold > 0, "Threshold must be > 0"
        assert self.block_size_threshold <= 1.0, "Threshold must be < 1"
        nframes = sorted([b[5] - b[4] for b in self.reflections['bbox']])
        cutoff = int(self.block_size_threshold*len(nframes))
        block_size = nframes[cutoff] * 2
        self.block_size = block_size
        self.block_size_units = 'frames'

  def compute_jobs(self):
    '''
    Compute the jobs

    '''
    from itertools import groupby
    from math import ceil
    groups = groupby(
      range(len(self.experiments)),
      lambda x: (self.experiments[x].imageset,
                 self.experiments[x].scan))
    self.jobs = JobList()
    for key, indices in groups:
      indices = list(indices)
      i0 = indices[0]
      i1 = indices[-1]+1
      expr = self.experiments[i0]
      scan = expr.scan
      imgs = expr.imageset
      array_range = (0, len(imgs))
      if scan is not None:
        assert len(imgs) == len(scan), "Invalid scan range"
        array_range = scan.get_array_range()
      if self.block_size is None:
        block_size_frames = array_range[1] - array_range[0]
      elif self.block_size_units == 'radians':
        phi0, dphi = scan.get_oscillation(deg=False)
        block_size_frames = int(ceil(self.block_size / dphi))
      elif self.block_size_units == 'degrees':
        phi0, dphi = scan.get_oscillation()
        block_size_frames = int(ceil(self.block_size / dphi))
      elif self.block_size_units == 'frames':
        block_size_frames = int(ceil(self.block_size))
      else:
        raise RuntimeError('Unknown block_size_units = %s' % block_size_units)
      self.jobs.add((i0, i1), array_range, block_size_frames)
    assert len(self.jobs) > 0, "Invalid number of jobs"

  def split_reflections(self):
    '''
    Split the reflections into partials or over job boundaries

    '''
    from logging import info

    # Optionally split the reflection table into partials, otherwise,
    # split over job boundaries
    if self.partials:
      num_full = len(self.reflections)
      self.reflections.split_partials()
      num_partial = len(self.reflections)
      assert num_partial >= num_full, "Invalid number of partials"
      if (num_partial > num_full):
        info(' Split %d reflections into %d partial reflections\n' % (
          num_full,
          num_partial))
    else:
      num_full = len(self.reflections)
      self.jobs.split(self.reflections)
      num_partial = len(self.reflections)
      assert num_partial >= num_full, "Invalid number of partials"
      if (num_partial > num_full):
        num_split = num_partial - num_full
        info(' Split %d reflections overlapping job boundaries\n' % num_split)

    # Compute the partiality
    self.reflections.compute_partiality(
      self.experiments,
      self.profile_model)

  def compute_processors(self):
    '''
    Compute the number of processors

    '''
    from libtbx.introspection import machine_memory_info
    from math import floor
    from dials.array_family import flex

    # Set the memory usage per processor
    if (self.mp_method == 'multiprocessing' and self.mp_nproc > 1):

      # Get the maximum shoebox memory
      max_memory = flex.max(self.jobs.shoebox_memory(
        self.reflections, self.flatten))

      # Compute percentage of max available. The function is not portable to
      # windows so need to add a check if the function fails. On windows no
      # warning will be printed
      memory_info = machine_memory_info()
      total_memory = memory_info.memory_total()
      if total_memory is not None:
        assert total_memory > 0, "Your system appears to have no memory!"
        limit_memory = total_memory * self.max_memory_usage
        njobs = int(floor(limit_memory / max_memory))
        if njobs < 1:
          raise RuntimeError('''
            No enough memory to run integration jobs. Possible solutions
            include increasing the percentage of memory allowed for shoeboxes or
            decreasing the block size.
              Total system memory: %g GB
              Limit shoebox memory: %g GB
              Max shoebox memory: %g GB
          ''' % (total_memory/1e9, limit_memory/1e9, max_memory/1e9))
        else:
          self.mp_nproc = min(self.mp_nproc, njobs)
          self.max_memory_usage = self.max_memory_usage / self.mp_nproc

  def summary(self):
    '''
    Get a summary of the processing

    '''
    from libtbx.table_utils import format as table

    # Compute the task table
    if self.experiments.all_stills():
      rows = [["#",
               "Group",
               "Frame From",
               "Frame To",
               "# Reflections"]]
      for i in range(len(self)):
        job = self.manager.job(i)
        group = job.index()
        f0, f1 = job.frames()
        n = self.manager.num_reflections(i)
        rows.append([str(i), str(group), str(f0), str(f1), str(n)])
    elif self.experiments.all_sweeps():
      rows = [["#",
               "Group",
               "Frame From",
               "Frame To",
               "Angle From",
               "Angle To",
               "# Reflections"]]
      for i in range(len(self)):
        job = self.manager.job(i)
        group = job.index()
        expr = job.expr()
        f0, f1 = job.frames()
        scan = self.experiments[expr[0]].scan
        p0 = scan.get_angle_from_array_index(f0)
        p1 = scan.get_angle_from_array_index(f1)
        n = self.manager.num_reflections(i)
        rows.append([str(i), str(group), str(f0), str(f1), str(p0), str(p1), str(n)])
    else:
      raise RuntimeError('Experiments must be all sweeps or all stills')

    # The job table
    task_table = table(rows, has_header=True, justify="right", prefix=" ")

    # The format string
    if self.block_size is None:
      block_size = "auto"
    else:
      block_size = str(self.block_size)
    fmt = (
      'Processing reflections in the following blocks of images:\n'
      '\n'
      ' block_size: %s %s\n'
      '\n'
      '%s\n'
    )
    return fmt % (block_size, self.block_size_units, task_table)

class ManagerRot(Manager):
  ''' Specialize the manager for oscillation data using the oscillation pre and
  post processors. '''

  def __init__(self, experiments, profile_model, reflections, params):
    ''' Initialise the pre-processor, post-processor and manager. '''

    # Ensure we have the correct type of data
    if not experiments.all_sweeps():
      raise RuntimeError('''
        An inappropriate processing algorithm may have been selected!
         Trying to perform rotation processing when not all experiments
         are indicated as rotation experiments.
      ''')

    # Initialise the manager
    super(ManagerRot, self).__init__(
      experiments, profile_model, reflections, params)


class ManagerStills(Manager):
  ''' Specialize the manager for stills data using the stills pre and post
  processors. '''

  def __init__(self, experiments, profile_model, reflections, params):
    ''' Initialise the pre-processor, post-processor and manager. '''

    # Ensure we have the correct type of data
    if not experiments.all_stills():
      raise RuntimeError('''
        An inappropriate processing algorithm may have been selected!
         Trying to perform stills processing when not all experiments
         are indicated as stills experiments.
      ''')

    # Initialise the manager
    super(ManagerStills, self).__init__(
      experiments, profile_model, reflections, params)


class Processor3D(Processor):
  ''' Top level processor for 3D processing. '''

  def __init__(self, experiments, profile_model, reflections, params):
    ''' Initialise the manager and the processor. '''

    # Set some parameters
    params.shoebox.partials = False
    params.shoebox.flatten = False

    # Create the processing manager
    manager = ManagerRot(experiments, profile_model, reflections, params)

    # Initialise the processor
    super(Processor3D, self).__init__(manager)


class ProcessorFlat3D(Processor):
  ''' Top level processor for flat 2D processing. '''

  def __init__(self, experiments, profile_model, reflections, params):
    ''' Initialise the manager and the processor. '''

    # Set some parameters
    params.shoebox.flatten = True
    params.shoebox.partials = False

    # Create the processing manager
    manager = ManagerRot(experiments, profile_model, reflections, params)

    # Initialise the processor
    super(ProcessorFlat3D, self).__init__(manager)


class Processor2D(Processor):
  ''' Top level processor for 2D processing. '''

  def __init__(self, experiments, profile_model, reflections, params):
    ''' Initialise the manager and the processor. '''

    # Set some parameters
    params.shoebox.partials = True

    # Create the processing manager
    manager = ManagerRot(experiments, profile_model, reflections, params)

    # Initialise the processor
    super(Processor2D, self).__init__(manager)


class ProcessorSingle2D(Processor):
  ''' Top level processor for still image processing. '''

  def __init__(self, experiments, profile_model, reflections, params):
    ''' Initialise the manager and the processor. '''

    # Set some of the parameters
    params.block.size = 1
    params.block.units = 'frames'
    params.shoebox.partials = True
    params.shoebox.flatten = False

    # Create the processing manager
    manager = ManagerRot(experiments, profile_model, reflections,  params)

    # Initialise the processor
    super(ProcessorSingle2D, self).__init__(manager)


class ProcessorStills(Processor):
  ''' Top level processor for still image processing. '''

  def __init__(self, experiments, profile_model, reflections, params):
    ''' Initialise the manager and the processor. '''

    # Set some parameters
    params.block.size = 1
    params.block.units = 'frames'
    params.shoebox.partials = False
    params.shoebox.flatten = False

    # Create the processing manager
    manager = ManagerStills(experiments, profile_model, reflections, params)

    # Initialise the processor
    super(ProcessorStills, self).__init__(manager)
