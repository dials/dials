#
# integrator.py
#
#  Copyright (C) 2013 Diamond Light Source
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
  ''' Generate the phil scope. '''
  import dials.extensions
  from dials.interfaces import BackgroundIface
  from dials.interfaces import IntensityIface
  from dials.interfaces import CentroidIface

  phil_scope = phil.parse('''

    integration {

      include scope dials.data.lookup.phil_scope
      include scope dials.data.multiprocessing.phil_scope

      reference {

        select = False
          .type = bool
          .help = "Do selection of reference profiles (if false use indexed"
                  "reflections)"

        truncate = 0.9,0.01
          .type = floats(size=2,value_min=0.0, value_max=1.0)
          .help = "Remove a fraction of lowest and highest intensity"
                  "reflections for selection."
      }

      block {
        size = auto
          .type = float
          .help = "The block size in rotation angle (degrees)."

        max_size = 5
          .type = float(value_min=0.0)
          .help = "The maximum size (in degrees) of an integration block"

        threshold = 0.95
          .type = float(value_min=0.0, value_max=1.0)
          .help = "For block size auto the block size is calculated by sorting"
                  "reflections by the number of frames they cover and then"
                  "selecting the block size to be 2*nframes[threshold] such"
                  "that 100*threshold % of reflections are guarenteed to be"
                  "fully contained in 1 block"

        max_mem_usage = 0.5
          .type = float(value_min=0.0,value_max=1.0)
          .help = "The maximum percentage of total physical memory to use for"
                  "allocating shoebox arrays."
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

        include scope dials.algorithms.integration.filtering.phil_scope
      }

      debug {
        save_shoeboxes = False
          .type = bool
          .help = "Save shoeboxes after each integration task."
      }
    }
  ''', process_includes=True)
  main_scope = phil_scope.get_without_substitution("integration")
  assert(len(main_scope) == 1)
  main_scope = main_scope[0]
  main_scope.adopt_scope(BackgroundIface.phil_scope())
  main_scope.adopt_scope(IntensityIface.phil_scope())
  main_scope.adopt_scope(CentroidIface.phil_scope())
  return phil_scope

# The integration phil scope
phil_scope = generate_phil_scope()


class Integrator(object):
  ''' Integrator interface class. '''

  def __init__(self, manager, nthreads=1, nproc=1, mp_method='multiprocessing'):
    ''' Initialise the integrator.

    The integrator requires a manager class implementing the IntegratorManager
    interface. This class executes all the workers in separate threads and
    accumulates the results to expose to the user.

    Params:
      manager The integration manager
      nproc The number of processors
      mp_method The multiprocessing method

    '''
    assert(nthreads > 0)
    assert(nproc > 0)
    self._manager = manager
    self._nthreads = nthreads
    self._nproc = nproc
    self._mp_method = mp_method
    if nproc > 1 and nthreads > 1 and mp_method == 'multiprocessing':
      raise RuntimeError( '''
        At the moment, running multiple jobs in parallel on a single machine
        with more than 1 thread per jobs is not implemented. Either set the
        number of threads per job to be 1 or set the number of jobs to be 1
        or set the multiprocessing method to operate on a cluster
      ''')

  def integrate(self):
    ''' Do all the integration tasks.

    Returns
      The integration results

    '''
    from time import time
    from libtbx import easy_mp
    from libtbx.table_utils import format as table
    from logging import info
    import omptbx

    start_time = time()
    if self._nthreads > omptbx.omp_get_num_procs():
      nthreads = omptbx.omp_get_num_procs()
    else:
      nthreads = self._nthreads
    omptbx.omp_set_num_threads(nthreads)
    job.nthreads = nthreads
    self._manager.initialize()
    num_proc = len(self._manager)
    if self._nproc > 0:
      num_proc = min(num_proc, self._nproc)
    info(' Using %s with %d parallel job(s) and %d thread(s) per job\n' % (
      self._mp_method, num_proc, nthreads))
    if num_proc > 1:
      def process_output(result):
        import logging
        for message in result[1]:
          logging.log(message.levelno, message.msg)
        self._manager.accumulate(result[0])
      def execute_task(task):
        from cStringIO import StringIO
        from dials.util import log
        import logging
        log.config_simple_cached()
        result = task()
        handlers = logging.getLogger().handlers
        assert(len(handlers) == 1)
        return result, handlers[0].messages()
      task_results = easy_mp.parallel_map(
        func=execute_task,
        iterable=list(self._manager.tasks()),
        processes=num_proc,
        callback=process_output,
        method=self._mp_method,
        preserve_order=True,
        preserve_exception_message=True)
      task_results, output = zip(*task_results)
    else:
      task_results = [task() for task in self._manager.tasks()]
      for result in task_results:
        self._manager.accumulate(result)
    self._manager.finalize()
    end_time = time()
    rows = [
      ["Read time"        , "%.2f seconds" % (self._manager.time().read)       ],
      ["Extract time"     , "%.2f seconds" % (self._manager.time().extract)    ],
      ["Pre-process time" , "%.2f seconds" % (self._manager.time().preprocess) ],
      ["Process time"     , "%.2f seconds" % (self._manager.time().process)    ],
      ["Post-process time", "%.2f seconds" % (self._manager.time().postprocess)],
      ["Total time"       , "%.2f seconds" % (self._manager.time().total)      ],
      ["User time"        , "%.2f seconds" % (end_time - start_time)           ],
    ]
    info('')
    info(table(rows, justify='right', prefix=' '))
    return self._manager.result()

def hist(data, width=80, symbol='#', prefix=''):
  ''' A utility function to print a histogram of reflections on frames. '''
  from collections import defaultdict, Counter
  from math import log10, floor
  assert(len(data) > 0)
  assert(width > 0)
  count = Counter(data)
  count = count.items()
  count.sort()
  frame, count = zip(*count)
  min_frame = min(frame)
  max_frame = max(frame)
  min_count = min(count)
  max_count = max(count)
  assert(max_count > 0)
  assert(min_count >= 0)
  if max_frame == 0:
    num_frame_zeros = 1
  else:
    num_frame_zeros = int(floor(log10(max_frame))) + 1
  num_count_zeros = int(floor(log10(max_count))) + 1
  assert(num_frame_zeros > 0)
  assert(num_count_zeros > 0)
  num_hist = width - (num_frame_zeros + num_count_zeros + 5) - len(prefix)
  assert(num_hist) > 0
  fmt = '%s%%-%dd [%%-%dd]: %%s' % (prefix, num_frame_zeros, num_count_zeros)
  scale = float(num_hist) / max_count
  return '\n'.join((
    fmt % (key, value, int(value * scale) * symbol)
      for key, value in zip(frame, count)))

def frame_hist(bbox, width=80, symbol='#', prefix=''):
  ''' A utility function to print a histogram of reflections on frames. '''
  return hist(
    [z for b in bbox for z in range(b[4], b[5])],
    width=width,
    symbol=symbol,
    prefix=prefix)

def nframes_hist(bbox, width=80, symbol='#', prefix=''):
  ''' A utility function to print a histogram of number of frames. '''
  return hist(
    [b[5] - b[4] for b in bbox],
    width=width,
    symbol=symbol,
    prefix=prefix)


class Result(object):
  ''' A class representing an integration result. '''

  def __init__(self, index, reflections):
    ''' Initialise the data. '''
    self.index = index
    self.reflections = reflections


class Task(object):
  ''' A class to perform an integration task. '''

  def __init__(self,
               index,
               experiments,
               profile_model,
               data,
               job,
               flatten=False,
               save_shoeboxes=False,
               max_mem_usage=0.5,
               reference_selector=None):
    ''' Initialise the task. '''
    self._index = index
    self._experiments = experiments
    self._profile_model = profile_model
    self._data = data
    self._job = job
    self._flatten = flatten
    self._save_shoeboxes = save_shoeboxes
    self._max_mem_usage = max_mem_usage
    self._mask = None
    self._reference_selector = reference_selector

  def __call__(self):
    ''' Do the integration. '''
    from dials.array_family import flex
    from time import time
    from dials.util.command_line import heading
    from libtbx.introspection import machine_memory_info
    from libtbx.table_utils import format as table
    from dials.algorithms.integration import statistics
    from logging import info

    # Get the start time
    start_time = time()

    # Set the global process ID
    job.index = self._index

    # Compute the number of shoebox bytes
    x0, x1, y0, y1, z0, z1 = self._data['bbox'].parts()
    assert(x1.all_gt(x0))
    assert(y1.all_gt(y0))
    assert(z1.all_gt(z0))
    xs = x1 - x0
    ys = y1 - y0
    zs = z1 - z0
    size = flex.sum(xs * ys * zs)
    nbytes = size * (4 + 4 + 4)

    # Compute percentage of max available. The function is not portable to
    # windows so need to add a check if the function fails. On windows no
    # warning will be printed
    memory_info = machine_memory_info()
    total_memory = memory_info.memory_total()
    if total_memory is not None:
      assert(total_memory > 0)
      assert(self._max_mem_usage >  0.0)
      assert(self._max_mem_usage <= 1.0)
      limit_memory = total_memory * self._max_mem_usage
      if nbytes > limit_memory:
        raise RuntimeError('''
        There was a problem allocating memory for shoeboxes. Possible solutions
        include increasing the percentage of memory allowed for shoeboxes or
        decreasing the block size.
          Total system memory: %g GB
          Limit shoebox memory: %g GB
          Requested shoebox memory: %g GB
        ''' % (total_memory/1e9, limit_memory/1e9, nbytes/1e9))

    # Print out some info
    EPS = 1e-7
    full_value = (0.997300203937 - EPS)
    fully_recorded = self._data['partiality'] > full_value
    num_partial = fully_recorded.count(False)
    num_full = fully_recorded.count(True)
    num_integrate = self._data.get_flags(self._data.flags.dont_integrate).count(False)
    num_reference = self._data.get_flags(self._data.flags.reference_spot).count(True)
    num_ice_ring = self._data.get_flags(self._data.flags.in_powder_ring).count(True)
    info("=" * 80)
    info("")
    info(heading("Beginning integration job %d" % self._index))
    info("")
    info(" Frames: %d -> %d" % self._job)
    info("")
    info(" Number of reflections")
    info("  Partial:     %d" % num_partial)
    info("  Full:        %d" % num_full)
    info("  In ice ring: %d" % num_ice_ring)
    info("  Integrate:   %d" % num_integrate)
    info("  Reference:   %d" % num_reference)
    info("  Total:       %d" % len(self._data))
    info("")
    info(" Memory usage")
    info("  Total system memory: %g GB" % (total_memory/1e9))
    info("  Limit shoebox memory: %g GB" % (limit_memory/1e9))
    info("  Requested shoebox memory: %g GB" % (nbytes/1e9))
    info("")

    # Get the sub imageset
    imagesets = self._experiments.imagesets()
    assert(len(imagesets) == 1)
    imageset = imagesets[0]
    frame00, frame01 = self._job
    try:
      frame10, frame11 = imageset.get_array_range()
    except Exception:
      frame10, frame11 = (0, len(imageset))
    assert(frame00 < frame01)
    assert(frame10 < frame11)
    assert(frame00 >= frame10)
    assert(frame01 <= frame11)
    index0 = frame00 - frame10
    index1 = index0 + (frame01 - frame00)
    assert(index0 < index1)
    assert(index0 >= 0)
    assert(index1 <= len(imageset))
    imageset = imageset[index0:index1]

    # Print a histogram of reflections on frames
    if frame01 - frame00 > 1:
      info(' The following histogram shows the number of reflections predicted')
      info(' to have all or part of their intensity on each frame.')
      info('')
      info(frame_hist(self._data['bbox'], prefix=' ', symbol='*'))
      info('')

    # Create the shoeboxes
    self._data["shoebox"] = flex.shoebox(
      self._data["panel"],
      self._data["bbox"],
      allocate=True)

    # Extract the shoeboxes
    read_time, extract_time = self._data.extract_shoeboxes(
      imageset,
      self._mask,
      job.nthreads)

    # Optionally flatten the shoeboxes
    if self._flatten:
      self._data['shoebox'].flatten()

    # Check if pixels are overloaded
    self._data.is_overloaded(self._experiments)

    # Process the data
    process_start_time = time()
    self._process()
    process_time = time() - process_start_time

    # Optionally save the shoeboxes
    if self._save_shoeboxes:
      filename = 'shoeboxes_%d.pickle' % self._index
      info('Saving shoeboxes to %s\n' % filename)
      self._data.as_pickle(filename)

    # Delete the shoeboxes
    del self._data["shoebox"]

    # Print some stats
    for stats in statistics.statistics(
        self._data,
        self._experiments,
        image_summary=False,
        resolution_summary=False,
        whole_summary=True):
      info(stats)

    # The total time
    total_time = time() - start_time

    # Print a table timing info for the task
    info('')
    info(table(
      [
        ["Read time"    , "%.2f seconds" % read_time    ],
        ["Extract time" , "%.2f seconds" % extract_time ],
        ["Process time" , "%.2f seconds" % process_time ],
        ["Total time"   , "%.2f seconds" % total_time   ],
      ],
      justify='right',
      prefix=' '))

    # Return the result
    result = Result(self._index, self._data)
    result.read_time = read_time
    result.extract_time = extract_time
    result.process_time = process_time
    result.total_time = total_time
    return result

  def _process(self):
    ''' Process the data. '''
    from logging import info
    self._data.compute_mask(self._experiments, self._profile_model)
    self._data.integrate(
      self._experiments,
      self._profile_model,
      reference_selector=self._reference_selector)
    info('')



class Manager(object):
  ''' An class to manage integration book-keeping '''

  class TimingInfo(object):
    ''' A class to contain timing info. '''
    def __init__(self):
      self.read = 0
      self.extract = 0
      self.preprocess = 0
      self.process = 0
      self.postprocess = 0
      self.total = 0

  def __init__(self,
               preprocess,
               postprocess,
               experiments,
               profile_model,
               reflections,
               block_size=libtbx.Auto,
               block_size_units='degrees',
               block_size_threshold=0.99,
               block_size_max=None,
               flatten=False,
               save_shoeboxes=False,
               max_mem_usage=0.5,
               reference_selector=None):
    ''' Initialise the manager. '''
    self._finalized = False
    self._flatten = flatten
    self._block_size = block_size
    self._block_size_units = block_size_units
    self._block_size_threshold = block_size_threshold
    self._block_size_max = block_size_max
    self._save_shoeboxes = save_shoeboxes
    self._max_mem_usage = max_mem_usage
    self._reference_selector = reference_selector
    self._preprocess = preprocess
    self._postprocess = postprocess
    self._experiments = experiments
    self._profile_model = profile_model
    self._reflections = reflections
    self._time = Manager.TimingInfo()

  def initialize(self):
    ''' Enter the context manager. '''
    from itertools import groupby
    from math import ceil
    self._compute_block_size()
    groups = groupby(
      range(len(self._experiments)),
      lambda x: (self._experiments[x].imageset,
                 self._experiments[x].scan))
    jobs = JobList()
    for key, indices in groups:
      indices = list(indices)
      i0 = indices[0]
      i1 = indices[-1]+1
      expr = self._experiments[i0]
      scan = expr.scan
      imgs = expr.imageset
      array_range = (0, len(imgs))
      if scan is not None:
        assert(len(imgs) == len(scan))
        array_range = scan.get_array_range()
      if self._block_size_units == 'degrees':
        phi0, dphi = scan.get_oscillation()
        block_size_frames = int(ceil(self._block_size / dphi))
      elif self._block_size_units == 'frames':
        block_size_frames = int(ceil(self._block_size))
      else:
        raise RuntimeError('Unknown block_size_units = %s' % block_size_units)
      jobs.add((i0, i1), array_range, block_size_frames)
    assert(len(jobs) > 0)
    self._preprocess(self._reflections, jobs)
    self._manager = ReflectionManager(jobs, self._reflections)
    self._time.preprocess = self._preprocess.time
    self._print_summary(self._block_size, self._block_size_units)

  def task(self, index):
    ''' Get a task. '''
    job = self._manager.job(index)
    frames = job.frames()
    expr_id = job.expr()
    assert(expr_id[1] > expr_id[0])
    assert(expr_id[0] >= 0)
    assert(expr_id[1] <= len(self._experiments))
    expr = self._experiments[expr_id[0]:expr_id[1]]
    prof = self._profile_model[expr_id[0]:expr_id[1]]
    data = self._manager.split(index)
    return Task(
      index=index,
      experiments=expr,
      profile_model=prof,
      data=data,
      job=frames,
      flatten=self._flatten,
      save_shoeboxes=self._save_shoeboxes,
      max_mem_usage=self._max_mem_usage,
      reference_selector=self._reference_selector)

  def tasks(self):
    ''' Iterate through the tasks. '''
    for i in range(len(self)):
      yield self.task(i)

  def accumulate(self, result):
    ''' Accumulate the results. '''
    self._manager.accumulate(result.index, result.reflections)
    self._time.read += result.read_time
    self._time.extract += result.extract_time
    self._time.process += result.process_time
    self._time.total += result.total_time

  def finalize(self):
    ''' Do the post-processing and finish. '''
    from dials.algorithms.integration import statistics
    from logging import info
    assert(self._manager.finished())
    self._postprocess(self._manager.data())
    self._time.postprocess = self._postprocess.time
    info('=' * 80)
    info('\n')
    for stats in statistics.statistics(
        self._manager.data(),
        self._experiments):
      info(stats)
    self._finalized = True

  def result(self):
    ''' Return the result. '''
    assert(self._finalized == True)
    return self._manager.data()

  def finished(self):
    ''' Return if all tasks have finished. '''
    return self._finalized and self._manager.finished()

  def __len__(self):
    ''' Return the number of tasks. '''
    return len(self._manager)

  def time(self):
    ''' Return the timing information. '''
    return self._time

  def _compute_block_size(self):
    ''' Compute the integration block size. '''
    from logging import info
    from math import ceil
    if self._block_size == libtbx.Auto:
      assert(self._block_size_threshold > 0)
      assert(self._block_size_threshold <= 1.0)
      self._reflections.compute_bbox(
        self._experiments,
        self._profile_model)
      nframes = sorted([b[5] - b[4] for b in self._reflections['bbox']])
      cutoff = int(self._block_size_threshold*len(nframes))
      block_size = nframes[cutoff] * 2
      block_size_max = 1
      for e in self._experiments:
        if e.scan:
          phi0, dphi = e.scan.get_oscillation()
          bm = int(ceil(self._block_size_max / dphi))
          if bm > block_size_max:
            block_size_max = bm
      if block_size > block_size_max:
        import bisect
        ind = bisect.bisect_right(nframes, block_size_max/2.0)
        perc = 100.0*(ind-1.0)/len(nframes)
        info('*' * 80)
        info('WARNING: setting block size to %d frames')
        info(' This may result in reflections being split across jobs')
        info(' automatic block size = %d frames' % block_size)
        info(' maximum block size = %d frames' % block_size_max)
        info(' %.2f%% of reflections guaranteed to be in a single block' % perc)
        info('*' * 80)
        block_size = block_size_max
      self._block_size = block_size
      self._block_size_units = 'frames'
    elif self._block_size_max is not None:
      if self._block_size_units == 'frames':
        assert(self._block_size < self._block_size_max)
      elif self._block_size_units == 'degrees':
        for e in self._experiments:
          phi0, dphi = e.scan.get_oscillation()
          bm = int(ceil(self._block_size_max / dphi))
          assert(self._block_size < self._block_size_max)
      else:
        raise RuntimeError('Unknown block_size_units = %s' % block_size_units)

  def _print_summary(self, block_size, block_size_units):
    ''' Print a summary of the integration stuff. '''
    from libtbx.table_utils import format as table
    from dials.util.command_line import heading
    from logging import info

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
      raise RuntimeError('Experiments must be all sweeps or all stills')
    task_table = table(rows, has_header=True, justify="right", prefix=" ")

    # Create summary format
    summary_format_str = (
      '%s\n'
      '\n'
      + heading('Beginning integration of the following experiments:\n') +
      '\n'
      ' Experiments: %d\n'
      ' Beams:       %d\n'
      ' Detectors:   %d\n'
      ' Goniometers: %d\n'
      ' Scans:       %d\n'
      ' Crystals:    %d\n'
      ' Imagesets:   %d\n'
      '\n'
      'Integrating reflections in the following blocks of images:\n'
      '\n'
      ' block_size: %d %s\n'
      '\n'
      '%s\n'
    )

    # Print the summary
    info(summary_format_str % (
      '=' * 80,
      len(self._experiments),
      len(self._experiments.beams()),
      len(self._experiments.detectors()),
      len(self._experiments.goniometers()),
      len(self._experiments.scans()),
      len(self._experiments.crystals()),
      len(self._experiments.imagesets()),
      block_size,
      block_size_units,
      task_table))


class PreProcessorRot(object):
  ''' A pre-processing class for oscillation data. '''

  def __init__(self,
               experiments,
               profile_model,
               min_zeta,
               powder_filter,
               partials,
               max_shoebox_overlap=1.0,
               **kwargs):
    ''' Initialise the pre-processor. '''
    self.experiments = experiments
    self.profile_model = profile_model
    self.min_zeta = min_zeta
    self.powder_filter = powder_filter
    self.partials = partials
    self.max_shoebox_overlap = max_shoebox_overlap
    self.time = 0

  def __call__(self, data, jobs):
    ''' Do some pre-processing. '''
    from dials.array_family import flex
    from dials.util.command_line import heading
    from scitbx.array_family import shared
    from logging import info
    from time import time
    info('=' * 80)
    info('')
    info(heading('Pre-processing reflections'))
    info('')
    st = time()

    # Compute some reflection properties
    data.compute_zeta_multi(self.experiments)
    data.compute_d(self.experiments)
    data.compute_bbox(self.experiments, self.profile_model)

    # Optionally split the reflection table into partials
    if self.partials:
      num_full = len(data)
      data.split_partials()
      num_partial = len(data)
      assert(num_partial >= num_full)
      if (num_partial > num_full):
        info(' Split %d reflections into %d partial reflections\n' % (
          num_full,
          num_partial))
    else:
      num_full = len(data)
      jobs.split(data)
      num_partial = len(data)
      assert(num_partial >= num_full)
      if (num_partial > num_full):
        num_split = num_partial - num_full
        info(' Split %d reflections overlapping job boundaries\n' % num_split)

    # Compute the partiality
    data.compute_partiality(self.experiments, self.profile_model)

    # Filter the reflections by zeta
    mask = flex.abs(data['zeta']) < self.min_zeta
    num_ignore = mask.count(True)
    data.set_flags(mask, data.flags.dont_integrate)

    # Filter the reflections by powder rings
    if self.powder_filter is not None:
      mask = self.powder_filter(data['d'])
      data.set_flags(mask, data.flags.in_powder_ring)

    # Print some output
    EPS = 1e-7
    full_value = (0.997300203937 - EPS)
    fully_recorded = data['partiality'] > full_value
    num_partial = fully_recorded.count(False)
    num_full = fully_recorded.count(True)
    num_integrate = data.get_flags(data.flags.dont_integrate).count(False)
    num_reference = data.get_flags(data.flags.reference_spot).count(True)
    num_ice_ring = data.get_flags(data.flags.in_powder_ring).count(True)

    # Find any overlaps
    overlaps = data.find_overlaps(self.experiments)
    num_overlap = 0
    for i in range(len(data)):
      if len(list(overlaps.adjacent_vertices(i))) > 0:
        num_overlap += 1
    perc_overlap = 100.0 * num_overlap / len(data)

    # Compute the fractional overlap
    # frac_overlap = data.compute_shoebox_overlap_fraction(overlaps)
    # mask = (frac_overlap > self.max_shoebox_overlap)
    # num_bad_overlap = mask.count(True)
    # data.set_flags(mask, data.flags.dont_integrate)


    self.time = time() - st
    info(' Number of reflections')
    info('  Partial:     %d' % num_partial)
    info('  Full:        %d' % num_full)
    info('  In ice ring: %d' % num_ice_ring)
    info('  Integrate:   %d' % num_integrate)
    info('  Reference:   %d' % num_reference)
    info('  Ignore:      %d' % num_ignore)
    info('  Total:       %d' % len(data))
    info('')
    info(' Filtered %d reflections by zeta = %0.3f' % (num_ignore, self.min_zeta))
    info('')
    info(' Shoebox overlaps (foreground and background)')
    info('  Reflections with 1 or more overlap: %d (%.2f%%)' % (num_overlap, perc_overlap))
    info('  Total number of overlaps:           %d' % overlaps.num_edges())
    info('')
    # info(' Filtered %d reflections by max_shoebox_overlap = %0.3f' %
    #      (num_bad_overlap, self.max_shoebox_overlap))
    info('')
    info(' Time taken: %.2f seconds' % self.time)
    info('')


class PreProcessorStills(object):
  ''' A pre-processing class for stills data. '''

  def __init__(self,
               experiments,
               profile_model,
               powder_filter,
               **kwargs):
    ''' Initialise the pre-processor. '''
    self.experiments = experiments
    self.profile_model = profile_model
    self.powder_filter = powder_filter
    self.time = 0

  def __call__(self, data, jobs):
    ''' Do some pre-processing. '''
    from dials.array_family import flex
    from dials.util.command_line import heading
    from logging import info
    from time import time
    info('=' * 80)
    info('')
    info(heading('Pre-processing reflections'))
    info('')
    st = time()

    # Compute some reflection properties
    data.compute_d(self.experiments)
    data.compute_bbox(self.experiments, self.profile_model)

    # Check the bounding boxes are all 1 frame in width
    z0, z1 = data['bbox'].parts()[4:6]
    assert((z1 - z0).all_eq(1))

    # Compute the partiality
    data.compute_partiality(self.experiments, self.profile_model)

    # Filter the reflections by powder rings
    if self.powder_filter is not None:
      mask = self.powder_filter(data['d'])
      data.set_flags(mask, data.flags.in_powder_ring)

    # Print out the pre-processing summary
    num_ignore = 0
    EPS = 1e-7
    full_value = (0.997300203937 - EPS)
    fully_recorded = data['partiality'] > full_value
    num_partial = fully_recorded.count(False)
    num_full = fully_recorded.count(True)
    num_integrate = data.get_flags(data.flags.dont_integrate).count(False)
    num_reference = data.get_flags(data.flags.reference_spot).count(True)
    num_ice_ring = data.get_flags(data.flags.in_powder_ring).count(True)
    self.time = time() - st
    info(' Number of reflections')
    info('  Partial:     %d' % num_partial)
    info('  Full:        %d' % num_full)
    info('  In ice ring: %d' % num_ice_ring)
    info('  Integrate:   %d' % num_integrate)
    info('  Reference:   %d' % num_reference)
    info('  Ignore:      %d' % num_ignore)
    info('  Total:       %d' % len(data))
    info('')
    info(' Time taken: %.2f seconds' % self.time)
    info('')


class PostProcessorRot(object):
  ''' A post-processing class for oscillation data. '''

  def __init__(self, experiments):
    ''' Initialise the post processor. '''
    self.experiments = experiments
    self.time = 0

  def __call__(self, data):
    ''' Do some post processing. '''
    from dials.util.command_line import heading
    from logging import info
    from time import time
    st = time()
    info('=' * 80)
    info('')
    info(heading('Post-processing reflections'))
    info('')

    # Compute the corrections
    data.compute_corrections(self.experiments)

    info('')
    num = data.get_flags(data.flags.integrated, all=False).count(True)
    self.time = time() - st
    info(' Integrated %d reflections' % num)
    info('')
    info(' Time taken: %.2f seconds' % self.time)
    info('')


class PostProcessorStills(object):
  ''' A post-processing class for stills data. '''

  def __init__(self, experiments):
    ''' Initialise the post processor. '''
    self.experiments = experiments
    self.time = 0

  def __call__(self, data):
    ''' Do some post processing. '''
    pass


class ManagerRot(Manager):
  ''' Specialize the manager for oscillation data using the oscillation pre and
  post processors. '''

  def __init__(self,
               experiments,
               profile_model,
               reflections,
               block_size=libtbx.Auto,
               block_size_units='degrees',
               block_size_threshold=0.99,
               block_size_max=None,
               min_zeta=0.05,
               powder_filter=None,
               flatten=False,
               partials=False,
               save_shoeboxes=False,
               max_mem_usage=0.5,
               reference_selector=None,
               max_shoebox_overlap=1.0,
               **kwargs):
    ''' Initialise the pre-processor, post-processor and manager. '''

    # Ensure we have the correct type of data
    if not experiments.all_sweeps():
      raise RuntimeError('''
        An inappropriate integration algorithm may have been selected!
         Trying to perform rotation integration when not all experiments
         are indicated as rotation experiments.
      ''')

    # Create the pre-processor
    preprocess = PreProcessorRot(
      experiments,
      profile_model,
      min_zeta=min_zeta,
      powder_filter=powder_filter,
      partials=partials,
      max_shoebox_overlap=max_shoebox_overlap)

    # Create the post-processor
    postprocess = PostProcessorRot(experiments)

    # Initialise the manager
    super(ManagerRot, self).__init__(
      preprocess,
      postprocess,
      experiments,
      profile_model,
      reflections,
      block_size=block_size,
      block_size_units=block_size_units,
      block_size_threshold=block_size_threshold,
      block_size_max=block_size_max,
      flatten=flatten,
      save_shoeboxes=save_shoeboxes,
      max_mem_usage=max_mem_usage,
      reference_selector=reference_selector)


class ManagerStills(Manager):
  ''' Specialize the manager for stills data using the stills pre and post
  processors. '''

  def __init__(self,
               experiments,
               profile_model,
               reflections,
               save_shoeboxes=False,
               max_mem_usage=0.5,
               powder_filter=None,
               reference_selector=None,
               **kwargs):
    ''' Initialise the pre-processor, post-processor and manager. '''

    # Ensure we have the correct type of data
    if not experiments.all_stills():
      raise RuntimeError('''
        An inappropriate integration algorithm may have been selected!
         Trying to perform stills integration when not all experiments
         are indicated as stills experiments.
      ''')

    # Create the pre-processor
    preprocess = PreProcessorStills(
      experiments,
      profile_model,
      powder_filter=powder_filter)

    # Create the post-processor
    postprocess = PostProcessorStills(experiments)

    # Initialise the manager
    super(ManagerStills, self).__init__(
      preprocess,
      postprocess,
      experiments,
      profile_model,
      reflections,
      block_size=1,
      block_size_units='frames',
      flatten=False,
      save_shoeboxes=save_shoeboxes,
      max_mem_usage=max_mem_usage,
      reference_selector=reference_selector)


class Integrator3D(Integrator):
  ''' Top level integrator for 3D integration. '''

  def __init__(self,
               experiments,
               profile_model,
               reflections,
               block_size=libtbx.Auto,
               block_size_threshold=0.99,
               block_size_max=None,
               min_zeta=0.05,
               powder_filter=None,
               nthreads=1,
               nproc=1,
               mp_method='multiprocessing',
               save_shoeboxes=False,
               max_mem_usage=0.5,
               reference_selector=None,
               max_shoebox_overlap=1.0,
               **kwargs):
    ''' Initialise the manager and the integrator. '''

    # Create the integration manager
    manager = ManagerRot(
      experiments,
      profile_model,
      reflections,
      block_size=block_size,
      block_size_threshold=block_size_threshold,
      block_size_max=block_size_max,
      min_zeta=min_zeta,
      max_shoebox_overlap=max_shoebox_overlap,
      powder_filter=powder_filter,
      save_shoeboxes=save_shoeboxes,
      max_mem_usage=max_mem_usage,
      reference_selector=reference_selector)

    # Initialise the integrator
    super(Integrator3D, self).__init__(manager, nthreads, nproc, mp_method)


class IntegratorFlat3D(Integrator):
  ''' Top level integrator for flat 2D integration. '''

  def __init__(self,
               experiments,
               profile_model,
               reflections,
               block_size=libtbx.Auto,
               block_size_threshold=0.99,
               block_size_max=None,
               min_zeta=0.05,
               max_shoebox_overlap=1.0,
               powder_filter=None,
               nthreads=1,
               nproc=1,
               mp_method='multiprocessing',
               save_shoeboxes=False,
               max_mem_usage=0.5,
               reference_selector=None,
               **kwargs):
    ''' Initialise the manager and the integrator. '''

    # Create the integration manager
    manager = ManagerRot(
      experiments,
      profile_model,
      reflections,
      block_size=block_size,
      block_size_threshold=block_size_threshold,
      block_size_max=block_size_max,
      min_zeta=min_zeta,
      max_shoebox_overlap=max_shoebox_overlap,
      powder_filter=powder_filter,
      flatten=True,
      save_shoeboxes=save_shoeboxes,
      max_mem_usage=max_mem_usage,
      reference_selector=reference_selector)

    # Initialise the integrator
    super(IntegratorFlat3D, self).__init__(manager, nthreads, nproc, mp_method)


class Integrator2D(Integrator):
  ''' Top level integrator for 2D integration. '''

  def __init__(self,
               experiments,
               profile_model,
               reflections,
               block_size=libtbx.Auto,
               block_size_threshold=0.99,
               block_size_max=None,
               min_zeta=0.05,
               max_shoebox_overlap=1.0,
               powder_filter=None,
               nthreads=1,
               nproc=1,
               mp_method='multiprocessing',
               save_shoeboxes=False,
               max_mem_usage=0.5,
               reference_selector=None,
               **kwargs):
    ''' Initialise the manager and the integrator. '''

    # Create the integration manager
    manager = ManagerRot(
      experiments,
      profile_model,
      reflections,
      block_size=block_size,
      block_size_threshold=block_size_threshold,
      block_size_max=block_size_max,
      min_zeta=min_zeta,
      max_shoebox_overlap=max_shoebox_overlap,
      powder_filter=powder_filter,
      partials=True,
      save_shoeboxes=save_shoeboxes,
      max_mem_usage=max_mem_usage,
      reference_selector=reference_selector)

    # Initialise the integrator
    super(Integrator2D, self).__init__(manager, nthreads, nproc, mp_method)


class IntegratorSingle2D(Integrator):
  ''' Top level integrator for still image integration. '''

  def __init__(self,
               experiments,
               profile_model,
               reflections,
               min_zeta=0.05,
               max_shoebox_overlap=1.0,
               powder_filter=None,
               nthreads=1,
               nproc=1,
               mp_method='multiprocessing',
               save_shoeboxes=False,
               max_mem_usage=0.5,
               reference_selector=None,
               **kwargs):
    ''' Initialise the manager and the integrator. '''

    # Create the integration manager
    manager = ManagerRot(
      experiments,
      profile_model,
      reflections,
      block_size=1,
      block_size_units='frames',
      min_zeta=min_zeta,
      max_shoebox_overlap=max_shoebox_overlap,
      powder_filter=powder_filter,
      partials=True,
      save_shoeboxes=save_shoeboxes,
      max_mem_usage=max_mem_usage,
      reference_selector=reference_selector)

    # Initialise the integrator
    super(IntegratorSingle2D, self).__init__(manager, nthreads, nproc, mp_method)


class IntegratorStills(Integrator):
  ''' Top level integrator for still image integration. '''

  def __init__(self,
               experiments,
               profile_model,
               reflections,
               powder_filter=None,
               nthreads=1,
               nproc=1,
               mp_method='multiprocessing',
               save_shoeboxes=False,
               max_mem_usage=0.5,
               reference_selector=None,
               **kwargs):
    ''' Initialise the manager and the integrator. '''

    # Create the integration manager
    manager = ManagerStills(
      experiments,
      profile_model,
      reflections,
      powder_filter=powder_filter,
      save_shoeboxes=save_shoeboxes,
      max_mem_usage=max_mem_usage,
      reference_selector=reference_selector)

    # Initialise the integrator
    super(IntegratorStills, self).__init__(manager, nthreads, nproc, mp_method)


class ReferenceSelector(object):
  ''' A class to select reflections for reference profiles. '''

  def __init__(self, truncate=(0.9,0.01)):
    ''' Set the parameters. '''
    self._truncate = tuple(truncate)

  def __call__(self, reflections):
    ''' Select the reference reflections. '''
    from dials.array_family import flex
    from logging import info
    from time import time
    st = time()
    minv = 100 * self._truncate[0]
    maxv = 100 * self._truncate[1]
    info('')
    info(' Selecting reference spots for profile fitting')
    info('  using %d reflections' % len(reflections))
    info('  removing %.1f%% of low intensity reflections' % minv)
    info('  removing %.1f%% of high intensity reflections' % maxv)

    # The reference flag
    reference_flag = reflections.flags.reference_spot

    # Unset the reference spot flags
    reflections.unset_flags(
      flex.size_t(range(len(reflections))),
      reference_flag)

    # Get the reflections integrated by summation
    integrated = reflections.get_flags(reflections.flags.integrated_sum)

    # Get the intensities and variances
    I = reflections['intensity.sum.value']

    # Get indices of those which are potential candidates
    indices = flex.size_t(range(len(I)))
    indices = indices.select(integrated)

    # Get the indices of reference spots
    indices = sorted(indices, key=lambda x: I[x])
    i0 = int(self._truncate[0] * len(indices))
    i1 = int(len(indices) - self._truncate[1] * len(indices))
    assert(i1 > i0)
    indices = flex.size_t(indices[i0:i1])

    # Get the min/max I/Sigma
    min_i = flex.min(I.select(indices))
    max_i = flex.max(I.select(indices))
    mean_i = flex.mean(I.select(indices))

    # Set the flags
    reflections.set_flags(indices, reference_flag)
    assert(reflections.get_flags(reference_flag).count(True) == len(indices))
    info('  selected %d reflections' % len(indices))
    info('  mean I: %g' % mean_i)
    info('  min  I: %g' % min_i)
    info('  max  I: %g' % max_i)
    info('  time taken: %g seconds' % (time() - st))

  @classmethod
  def from_params(Class, params):
    ''' Create the reference selector. '''
    if params.select:
      return Class(params.truncate)
    return None


class IntegratorFactory(object):
  ''' A factory for creating integrators. '''

  @staticmethod
  def create(params, experiments, profile_model, reflections):
    ''' Create the integrator from the input configuration. '''
    from dials.algorithms.integration.filtering import MultiPowderRingFilter
    from dials.interfaces import IntensityIface
    from dials.interfaces import BackgroundIface
    from dials.interfaces import CentroidIface
    from dials.array_family import flex
    from libtbx.utils import Abort

    # Check the input
    assert(len(experiments) == len(profile_model))

    # Check each experiment has an imageset
    for exp in experiments:
      if exp.imageset is None:
        raise Abort('''
          One or more experiment does not contain an imageset. Access to the
          image data is crucial for integration.
        ''')

    # Initialise the strategy classes
    BackgroundAlgorithm = BackgroundIface.extension(
      params.integration.background.algorithm)
    IntensityAlgorithm = IntensityIface.extension(
      params.integration.intensity.algorithm)
    CentroidAlgorithm = CentroidIface.extension(
      params.integration.centroid.algorithm)

    # Set the algorithms in the reflection table
    flex.reflection_table._background_algorithm = \
      flex.strategy(BackgroundAlgorithm, params)
    flex.reflection_table._intensity_algorithm = \
      flex.strategy(IntensityAlgorithm, params)
    flex.reflection_table._centroid_algorithm = \
      flex.strategy(CentroidAlgorithm, params)

    # Get the integrator class
    IntegratorClass = IntegratorFactory.select_integrator(
      IntensityAlgorithm.type(params, experiments))

    # Create the powder filter
    powder_filter = MultiPowderRingFilter.from_params(
      params.integration.filter)

    # Create the reference selector
    reference_selector = ReferenceSelector.from_params(
      params.integration.reference)

    # Compute the maximum memory usage per job
    max_mem_usage = params.integration.block.max_mem_usage
    if params.integration.mp.method == 'multiprocessing':
      max_mem_usage = max_mem_usage / params.integration.mp.nproc

    # Return an instantiation of the class
    return IntegratorClass(
      experiments=experiments,
      profile_model=profile_model,
      reflections=reflections,
      block_size=params.integration.block.size,
      block_size_threshold=params.integration.block.threshold,
      block_size_max=params.integration.block.max_size,
      min_zeta=params.integration.filter.min_zeta,
      max_shoebox_overlap=params.integration.filter.max_shoebox_overlap,
      powder_filter=powder_filter,
      nthreads=params.integration.mp.nthreads,
      nproc=params.integration.mp.nproc,
      mp_method=params.integration.mp.method,
      save_shoeboxes=params.integration.debug.save_shoeboxes,
      max_mem_usage=max_mem_usage,
      reference_selector=reference_selector)

  @staticmethod
  def select_integrator(integrator_type):
    ''' Select the integrator. '''
    if integrator_type == '3d':
      IntegratorClass = Integrator3D
    elif integrator_type == 'flat3d':
      IntegratorClass = IntegratorFlat3D
    elif integrator_type == '2d':
      IntegratorClass = Integrator2D
    elif integrator_type == 'single2d':
      IntegratorClass = IntegratorSingle2D
    elif integrator_type == 'stills':
      IntegratorClass = IntegratorStills
    else:
      raise RuntimeError("Unknown integration type")
    return IntegratorClass
