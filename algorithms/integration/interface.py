from __future__ import division
from dials.algorithms.integration import IntegrationTask3DExecutor
import boost.python
from iotbx import phil
import abc

# The integrator phil scope
phil_scope = phil.parse('''

  mp {
    method = *multiprocessing sge lsf pbs
      .type = choice
      .help = "The multiprocessing method to use"

    max_procs = 1
      .type = int(value_min=1)
      .help = "The number of processes to use."

    max_tasks = 1
      .type = int(value_min=1)
      .help = "The number of tasks to split the integration into. This number"
              "is a multiple of the maximum number of processes. This"
              "parameter controls how much the integration jobs (as specified"
              "by block.size are split up. Splitting up the jobs into fewer"
              "tasks makes the process more work efficient in terms of the"
              "number of images that need to be read; however, splitting the"
              "jobs into more tasks makes better use of the parallelism where"
              "some jobs take longer than others."
  }

  block {
    size = 10
    .type = float
    .help = "The block size in rotation angle (degrees)."
  }

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


class IntegrationTaskExecutor3DAux(boost.python.injector,
                                   IntegrationTask3DExecutor):
  ''' A class to add aditional methods to the executor class '''

  def execute(self, imageset, mask=None):
    ''' Passing in an imageset process all the images. '''
    from dials.model.data import Image
    import sys
    detector = imageset.get_detector()
    frame0 = self.frame0()
    frame1 = self.frame1()
    if mask is None:
      image = imageset[0]
      if not isinstance(image, tuple):
        image = (image,)
      mask = []
      for i in range(len(image)):
        tr = detector[i].get_trusted_range()
        mask.append(image[i].as_double() > tr[0])
      mask = tuple(mask)
    sys.stdout.write("Reading images: ")
    for image in imageset:
      if not isinstance(image, tuple):
        image = (image,)
      self.next(Image(image, mask))
      sys.stdout.write(".")
      sys.stdout.flush()
      del image
    sys.stdout.write("\n")
    sys.stdout.flush()
    assert(self.finished())


class IntegrationTask3D(IntegrationTask):
  ''' A class to perform a 3D integration task. '''

  def __init__(self, index, experiments, spec):
    ''' Initialise the task. '''
    self._index = index
    self._experiments = experiments
    self._spec = spec
    self._mask = None

  def __call__(self):
    ''' Do the integration. '''
    from dials.array_family import flex
    import sys
    print "=" * 80
    print ""
    print "Integrating task %d" % self._index

    class Process(object):

      def __init__(self):
        pass

      def __call__(self, reflections):
        print reflections
        return reflections

    process = Process()

    imageset = self._experiments[0].imageset

    executor = IntegrationTask3DExecutor(self._spec, process)
    executor.execute(imageset)

    return IntegrationResult(self._index, executor.data())


class IntegrationManager3D(IntegrationManager):
  ''' An class to manage 3D integration. book-keeping '''

  def __init__(self, experiments, reflections, params):
    ''' Initialise the manager. '''
    from dials.algorithms.integration import IntegrationManager3DExecutor
    from dials.array_family import flex
    imagesets = experiments.imagesets()
    detectors = experiments.detectors()
    scans = experiments.scans()
    assert(len(imagesets) == 1)
    assert(len(scans) == 1)
    assert(len(detectors) == 1)
    imageset = imagesets[0]
    scan = scans[0]
    detector = detectors[0]
    assert(len(imageset) == len(scan))
    self._experiments = experiments
    self._reflections = reflections
    phi0, dphi = scan.get_oscillation()
    block_size = params.block.size / dphi
    self._data = IntegrationManager3DExecutor(
      self._reflections,
      scan.get_array_range(),
      block_size,
      params.mp.max_procs * params.mp.max_tasks,
      len(detector))

  def task(self, index):
    ''' Get a task. '''
    return IntegrationTask3D(
      index,
      self._experiments,
      self._data.split(index))

  def tasks(self):
    ''' Iterate through the tasks. '''
    for i in range(len(self)):
      yield self.task(i)

  def accumulate(self, result):
    ''' Accumulate the results. '''
    self._data.accumulate(result.index, result.reflections)

  def result(self):
    ''' Return the result. '''
    return self._data.data()

  def finished(self):
    ''' Return if all tasks have finished. '''
    return self._data.finished()

  def __len__(self):
    ''' Return the number of tasks. '''
    return len(self._data)


# class IntegrationTask3D(IntegrationTask):
#   ''' A class to perform a 3D integration task. '''

#   def __init__(self, index, task, experiments, reflections):
#     ''' Initialise the task. '''
#     self._index = index
#     self._task = task
#     self._experiments = experiments
#     self._reflections = reflections

#   def __call__(self):
#     ''' Do the integration. '''
#     from dials.array_family import flex
#     import sys
#     import operator
#     from time import time
#     print "=" * 80
#     print ""
#     print "Integrating task %d" % self._index
#     self._reflections["shoebox"] = flex.shoebox(
#       self._reflections["panel"],
#       self._reflections["bbox"])
#     st = time()
#     self._reflections["shoebox"].allocate()
#     print time() - st
#     frame0, frame1 = self._task
#     imageset = self._experiments[0].imageset
#     st = time()
#     self._reflections.fill_shoeboxes(imageset[frame0:frame1])
#     print time() - st
#     del self._reflections["shoebox"]

#     return IntegrationResult(self._index, self._reflections)


# class IntegrationManager3D(IntegrationManager):
#   ''' An class to manage 3D integration. book-keeping '''

#   def __init__(self, experiments, reflections, params):
#     ''' Initialise the manager. '''
#     from dials.algorithms.integration import IntegrationManagerData3D
#     from dials.array_family import flex
#     block_size = params.block.size
#     imagesets = experiments.imagesets()
#     scans = experiments.scans()
#     assert(len(imagesets) == 1)
#     assert(len(scans) == 1)
#     imageset = imagesets[0]
#     scan = scans[0]
#     self._experiments = experiments
#     self._reflections = reflections
#     reflections.compute_zeta_multi(experiments)
#     reflections.compute_d(experiments)
#     self._data = IntegrationManagerData3D(
#       self._reflections,
#       scan.get_oscillation(deg=True),
#       scan.get_array_range(),
#       block_size)
#     self._print_summary(block_size)

#   def task(self, index):
#     ''' Get a task. '''
#     return IntegrationTask3D(
#       index,
#       self._data.block(index),
#       self._experiments,
#       self._data.split(index))

#   def tasks(self):
#     ''' Iterate through the tasks. '''
#     for i in range(len(self)):
#       yield self.task(i)

#   def accumulate(self, result):
#     ''' Accumulate the results. '''
#     self._data.accumulate(result.index, result.reflections)

#   def result(self):
#     ''' Return the result. '''
#     return self._data.data()

#   def finished(self):
#     ''' Return if all tasks have finished. '''
#     return self._data.finished()

#   def __len__(self):
#     ''' Return the number of tasks. '''
#     return len(self._data)

#   def _print_summary(self, block_size):
#     ''' Print a summary of the integration stuff. '''
#     from math import floor, log10
#     from libtbx import table_utils

#     # Create a table of integration tasks
#     rows = [
#       ["Frame From", "Frame To",
#        "Angle From", "Angle To",
#        "# Process",  "# Include"]
#     ]
#     scans = self._experiments.scans()
#     assert(len(scans) == 1)
#     for i in range(len(self)):
#       f0, f1 = self._data.block(i)
#       p0 = scans[0].get_angle_from_array_index(f0)
#       p1 = scans[0].get_angle_from_array_index(f1)
#       n0 = len(self._data.to_process(i))
#       n1 = len(self._data.to_include(i))
#       rows.append([str(f0), str(f1), str(p0), str(p1), str(n0), str(n1)])
#     task_table = table_utils.format(
#       rows,
#       has_header=True,
#       justify="right",
#       prefix=" ")

#     # Create summary format
#     summary_format_str = (
#       '%s\n'
#       '\n'
#       'Beginning integration of the following experiments:\n'
#       '\n'
#       ' Experiments: %d\n'
#       ' Beams:       %d\n'
#       ' Detectors:   %d\n'
#       ' Goniometers: %d\n'
#       ' Scans:       %d\n'
#       ' Crystals:    %d\n'
#       '\n'
#       'Integrating reflections in the following blocks of images:\n'
#       '\n'
#       ' block_size: %d degrees\n'
#       '\n'
#       '%s\n'
#       '\n'
#       ' %d reflections overlapping blocks removed from integration\n'
#       '\n'
#       ' If you see poor performance, it may be because DIALS is using too\n'
#       ' much memory. This could be because you have specified too many\n'
#       ' processes on the same machine, in which case try setting the number\n'
#       ' of processes to a smaller number. It could also be because the block\n'
#       ' size is too large, in which case try setting the block size to a\n'
#       ' smaller value. Reflections not fully recorded in a block are not\n'
#       ' integrated so be sure to check that in reducing the block size you\n'
#       ' are not throwing away too many reflections.\n'
#     )

#     # Print the summary
#     print summary_format_str % (
#       '=' * 80,
#       len(self._experiments),
#       len(self._experiments.beams()),
#       len(self._experiments.detectors()),
#       len(self._experiments.goniometers()),
#       len(self._experiments.scans()),
#       len(self._experiments.crystals()),
#       block_size,
#       task_table,
#       len(self._data.to_not_process()))


class Integrator3D(Integrator):
  ''' Top level integrator for 3D integration. '''

  def __init__(self, experiments, reflections, params):
    ''' Initialise the manager and the integrator. '''
    manager = IntegrationManager3D(experiments, reflections, params)
    super(Integrator3D, self).__init__(manager, params)
