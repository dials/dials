from __future__ import division
from iotbx import phil
import abc

# The integrator phil scope
phil_scope = phil.parse('''

  mp {
    method = *multiprocessing sge lsf pbs
      .type = choice
      .help = "The multiprocessing method to use"

    max_procs = 1
      .type = int(value_min=0)
      .help = "The number of processes to use"
  }

  block_size = 10
    .type = float
    .help = "The block size in rotation angle (degrees)."

  filter {

    min_zeta = 0.05
      .help = "Filter the reflections by the value of zeta. A value of less"
              "than or equal to zero indicates that this will not be used. A"
              "positive value is used as the minimum permissable value."
      .type = float

    ice_rings {
      filter = False
        .type = bool
      unit_cell = 4.498,4.498,7.338,90,90,120
        .type = unit_cell
        .help = "The unit cell to generate d_spacings for ice rings."
      space_group = 194
        .type = space_group
        .help = "The space group used to generate d_spacings for ice rings."
      d_min = 0
        .type = int(value_min=0)
        .help = "The minimum resolution to filter ice rings"
      width = 0.06
        .type = float(value_min=0.0)
        .help = "The width of an ice ring (in d-spacing)."
    }
  }
''')

class IntegrationResult(object):
  ''' A class representing an integration result. '''

  def __init__(self, index, reflections):
    ''' Initialise the data. '''
    self.index = index
    self.reflections = reflections


class IntegrationTask(object):
  ''' The integration task interface. '''

  __metaclass__ = abc.ABCMeta

  @abc.abstractmethod
  def __call__(self):
    pass


class IntegrationManager(object):
  ''' The integration manager interface. '''

  __metaclass__ = abc.ABCMeta

  @abc.abstractmethod
  def task(self, index):
    pass

  @abc.abstractmethod
  def tasks(self):
    pass

  @abc.abstractmethod
  def accumulate(self, result):
    pass

  @abc.abstractmethod
  def finished(self):
    pass

  @abc.abstractmethod
  def __len__(self):
    pass


class Integrator(object):
  ''' Integrator interface class. '''

  def __init__(self, manager, params):
    ''' Initialise the integrator.

    The integrator requires a manager class implementing the IntegratorManager
    interface. This class executes all the workers in separate threads and
    accumulates the results to expose to the user.

    Params:
      manager The integration manager
      params The parameters

    '''
    self._params = params
    self._manager = manager

  def integrate(self):
    ''' Do all the integration tasks.

    Returns
      The integration results

    '''
    from time import time
    from libtbx import easy_mp
    start_time = time()
    num_proc = len(self._manager)
    if self._params.mp.max_procs > 0:
      num_proc = min(num_proc, self._params.mp.max_procs)
    if num_proc > 1:
      def print_output(result):
        print result[1]
      def execute_task(task):
        from cStringIO import StringIO
        import sys
        sys.stdout = StringIO()
        result = task()
        output = sys.stdout.getvalue()
        return result, output
      task_results = easy_mp.parallel_map(
        func=execute_task,
        iterable=list(self._manager.tasks()),
        processes=num_proc,
        callback=print_output,
        method=self._params.mp.method,
        preserve_order=True)
      task_results, output = zip(*task_results)
    else:
      task_results = [task() for task in self._manager.tasks()]
    for result in task_results:
      self._manager.accumulate(result)
    end_time = time()
    print "Time taken: %.2f seconds" % (end_time - start_time)
    return self._manager.result()


class IntegrationTask3D(IntegrationTask):

  def __init__(self, index, experiments, task):
    self._index = index
    self._experiments = experiments
    self._task = task

  def __call__(self):

    print "=" * 80
    print ""
    print "Integrating task %d" % self._index
    imageset = self._experiments[0].imageset
    for i in range(*self._task):
      print "Reading Image % i" % i
      imageset[i]

    return IntegrationResult(self._index, None)


class IntegrationManager3D(IntegrationManager):

  def __init__(self, experiments, reflections, params):
    from dials.array_family import flex
    from dials.algorithms.integration import Preprocessor
    from dials.algorithms.integration import MultiPowderRingFilter
    from dials.algorithms.integration import PowderRingFilter
    block_size = params.block_size
    imagesets = experiments.imagesets()
    scans = experiments.scans()
    assert(len(imagesets) == 1)
    assert(len(scans) == 1)
    imageset = imagesets[0]
    scan = scans[0]
    self._experiments = experiments
    self._reflections = reflections
    self._tasks = self._compute_tasks(
      scan.get_oscillation(deg=True),
      scan.get_array_range(),
      block_size)
    self._finished = flex.bool(len(self._tasks), False)
    reflections.compute_zeta_multi(experiments)
    reflections.compute_d(experiments)
    self._print_summary(block_size)

  def task(self, index):
    assert(0 <= index < len(self))
    return IntegrationTask3D(
      index,
      self._experiments,
      self._tasks[index])

  def tasks(self):
    for i in range(len(self)):
      yield self.task(i)

  def accumulate(self, result):
    assert(0 <= result.index < len(self))
    self._finished[result.index] = True

  def result(self):
    assert(self.finished())
    return self._reflections

  def finished(self):
    return self._finished.all_eq(True)

  def __len__(self):
    return len(self._finished)

  def _compute_tasks(self, oscillation, array_range, block_size):
    ''' Compute the number of blocks. '''
    from math import ceil
    phi0, dphi = oscillation
    frame0, frame1 = array_range
    nframes = frame1 - frame0
    half_block_size = block_size / 2.0
    assert(nframes > 0)
    assert(half_block_size >= dphi)
    half_block_length = float(half_block_size) / dphi
    nblocks = int(ceil(nframes / half_block_length))
    assert(nblocks <= nframes)
    half_block_length = int(ceil(nframes / nblocks))
    blocks = [frame0]
    for i in range(nblocks):
      frame = frame0 + (i + 1) * half_block_length
      if frame > frame0+nframes:
        frame = frame0+nframes
      blocks.append(frame)
      if frame == frame0+nframes:
        break
    assert(all(b > a for a, b in zip(blocks, blocks[1:])))
    tasks = []
    for i in range(len(blocks)-2):
      i0 = blocks[i]
      i1 = blocks[i+2]
      tasks.append((i0, i1))
    return tasks

  def _print_summary(self, block_size):
    from math import floor, log10
    scans = self._experiments.scans()
    assert(len(scans) == 1)
    scan = scans[0]
    npad1 = int(floor(log10(max([max(t) for t in self._tasks])))) + 1
    npad2 = int(floor(log10(max(scan.get_oscillation_range())))) + 3
    format_string = ' %%%dd -> %%%dd (%%%d.2f -> %%%d.2f degrees)'
    format_string = format_string % (npad1, npad1, npad2, npad2)
    tasks_string = []
    for i in range(len(self._tasks)):
      f0, f1 = self._tasks[i]
      p0 = scan.get_angle_from_array_index(f0)
      p1 = scan.get_angle_from_array_index(f1)
      tasks_string.append(format_string % (f0, f1, p0, p1))
    tasks_string = '\n'.join(tasks_string)
    summary_format_str = (
      '%s\n'
      '\n'
      'Beginning integration of the following experiments:\n'
      '\n'
      ' Experiments: %d\n'
      ' Beams:       %d\n'
      ' Detectors:   %d\n'
      ' Goniometers: %d\n'
      ' Scans:       %d\n'
      ' Crystals:    %d\n'
      '\n'
      'Integrating reflections in the following blocks of images:\n'
      '\n'
      ' block_size: %d degrees\n'
      '\n'
      '%s\n'
      '\n'
    )
    print summary_format_str % (
      '=' * 80,
      len(self._experiments),
      len(self._experiments.beams()),
      len(self._experiments.detectors()),
      len(self._experiments.goniometers()),
      len(self._experiments.scans()),
      len(self._experiments.crystals()),
      block_size,
      tasks_string)


class Integrator3D(Integrator):
  ''' Top level integrator for 3D integration. '''

  def __init__(self, experiments, reflections, params):
    ''' Initialise the manager and the integrator. '''
    manager = IntegrationManager3D(experiments, reflections, params)
    super(Integrator3D, self).__init__(manager, params)
