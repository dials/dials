from __future__ import division
import abc

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

  def __init__(self, manager, mp_method='multiprocessing', max_procs=None):
    ''' Initialise the integrator.

    The integrator requires a manager class implementing the IntegratorManager
    interface. This class executes all the workers in separate threads and
    accumulates the results to expose to the user.

    Params:
      manager The integration manager
      mp_method The multi processing method
      max_procs The maximum number of processes

    '''
    self._mp_method = mp_method
    self._manager = manager
    self._max_procs = max_procs

  def integrate(self):
    ''' Do all the integration tasks.

    Returns
      The integration results

    '''
    from libtbx import easy_mp
    num_proc = len(self._manager)
    if self._max_procs is not None:
      num_proc = min(num_proc, self._max_procs)
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
        method=self._mp_method,
        preserve_order=True)
      task_results, output = zip(*task_results)
    else:
      task_results = [task() for task in self._manager.tasks()]
    for result in task_results:
      self._manager.accumulate(result)
    return self._manager.result()


# class ReflectionPreprocessor(object):

#   def __init__(self, reflections):
#     pass

#   def summary(self):
#     format_string = (
#       'Preprocessing reflections:\n'
#       '\n'
#       ' Number of reflections:                  %d\n'
#       ' Number of strong reflections:           %d\n'
#       ' Number filtered with zeta lt %.2f:      %d (%.2f%%)\n'
#       ' Number of reflection to integrate:      %d (%.2f%%)\n'
#       ' Number overlapping (background):        %d (%.2f%%)\n'
#       ' Number overlapping (foreground):        %d (%.2f%%)\n'
#       ' Number recorded on ice rings:           %d (%.2f%%)\n'
#       ' Number clipped at task boundaries:      %d (%.2f%%)\n'
#       ' Number strong overlapping (background): %d (%.2f%%)\n'
#       ' Number strong overlapping (foreground): %d (%.2f%%)\n'
#       ' Number strong recorded on ice rings:    %d (%.2f%%)\n'
#       ' Number strong for reference creation:   %d (%.2f%%)\n'
#       ' Smallest measurement box size (x):      %d\n'
#       ' Smallest measurement box size (y):      %d\n'
#       ' Smallest measurement box size (z):      %d\n'
#       ' Largest measurement box size (x):       %d\n'
#       ' Largest measurement box size (y):       %d\n'
#       ' Largest measurement box size (z):       %d\n'
#     )
#     return format_string % (0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
#                             0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)


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

    return IntegrationResult(self._index, None)


class IntegrationManager3D(IntegrationManager):

  def __init__(self, experiments, reflections, num_tasks=1, max_overlap=0):
    from dials.array_family import flex
    from dials.algorithms.integration import Preprocessor
    from dials.algorithms.integration import MultiPowderRingFilter
    from dials.algorithms.integration import PowderRingFilter
    imagesets = experiments.imagesets()
    scans = experiments.scans()
    assert(len(imagesets) == 1)
    assert(len(scans) == 1)
    imageset = imagesets[0]
    scan = scans[0]
    self._experiments = experiments
    self._reflections = reflections
    self._finished = flex.bool(num_tasks, False)
    self._tasks = self._compute_tasks(
      scan.get_array_range(),
      num_tasks,
      max_overlap)
    reflections.compute_zeta_multi(experiments)
    reflections.compute_d(experiments)
    filter = MultiPowderRingFilter()
    from cctbx.uctbx import unit_cell
    from cctbx.sgtbx import space_group_info
    # filter.add(PowderRingFilter(
    #   unit_cell((4.498,4.498,7.338,90,90,120)),
    #   space_group_info(194).group(),
    #   1,
    #   0.01))
    min_zeta = 0.05

    self._preprocessing = Preprocessor(reflections, filter, min_zeta)
    self._print_summary(num_tasks, max_overlap)

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

  def _compute_tasks(self, array_range, ntasks, max_overlap):
    ''' Compute the number of tasks. '''
    from math import ceil
    frame0, frame1 = array_range
    nframes = frame1 - frame0
    assert(nframes > 0)
    assert(ntasks <= nframes)
    task_length = int(ceil(nframes / ntasks))
    assert(task_length > 0)
    indices = [frame0]
    for i in range(ntasks):
      frame = frame0 + (i + 1) * task_length
      if frame > frame0+nframes:
        frame = frame0+nframes
      indices.append(frame)
      if frame == frame0+nframes:
        break
    assert(all(b > a for a, b in zip(indices, indices[1:])))
    tasks = []
    for i in range(len(indices)-1):
      i0 = indices[i]
      i1 = indices[i+1]
      i0 = max(frame0, i0-max_overlap)
      i1 = min(frame1, i1+max_overlap)
      tasks.append((i0, i1))
    for i in range(len(tasks)-1):
      i00, i01 = tasks[i]
      i10, i11 = tasks[i+1]
      assert(i00 < i10 and i01 < i11)
    return tasks

  def _print_summary(self, num_tasks, max_overlap):
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
      'Computing tasks based on the following parameters:\n'
      '\n'
      ' num_tasks:   %d\n'
      ' max_overlap: %d\n'
      '\n'
      'Integrating reflections in the following blocks of images:\n'
      '\n'
      '%s\n'
      '\n'
      '%s\n'
    )

    print summary_format_str % (
      '=' * 80,
      len(self._experiments),
      len(self._experiments.beams()),
      len(self._experiments.detectors()),
      len(self._experiments.goniometers()),
      len(self._experiments.scans()),
      len(self._experiments.crystals()),
      num_tasks,
      max_overlap,
      tasks_string,
      self._preprocessing.summary())


class Integrator3D(Integrator):
  ''' Top level integrator for 3D integration. '''

  def __init__(self,
               experiments,
               reflections,
               num_tasks=1,
               max_overlap=0,
               mp_method='multiprocessing',
               max_procs=None):

    # The 3D Integration manager
    manager = IntegrationManager3D(
      experiments=experiments,
      reflections=reflections,
      num_tasks=num_tasks,
      max_overlap=max_overlap)

    # Initialise the base integrator
    super(Integrator3D, self).__init__(
      manager=manager,
      mp_method=mp_method,
      max_procs=max_procs)
