from __future__ import division
from dials.algorithms.integration import IntegrationTask3DExecutor
from dials.algorithms.integration import IntegrationTask3DMultiExecutorBase
import boost.python
from iotbx import phil
import abc

def generate_phil_scope():
  ''' Generate the phil scope. '''
  from dials.interfaces import BackgroundIface
  from dials.interfaces import IntensityIface
  from dials.interfaces import CentroidIface

  phil_scope = phil.parse('''

    integration {

      include scope dials.data.lookup.phil_scope

      mp {
        method = *multiprocessing sge lsf pbs
          .type = choice
          .help = "The multiprocessing method to use"

        max_procs = 1
          .type = int(value_min=1)
          .help = "The number of processes to use."
      }

      block {
        size = 10
          .type = float
          .help = "The block size in rotation angle (degrees)."
      }

      shoebox {
        n_sigma = 3
          .help = "The number of standard deviations of the beam divergence and the"
                  "mosaicity to use for the bounding box size."
          .type = float

        sigma_b = None
          .help = "The E.S.D. of the beam divergence"
          .type = float

        sigma_m = None
          .help = "The E.S.D. of the reflecting range"
          .type = float
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

  def __init__(self, manager, max_procs=1, mp_method='multiprocessing'):
    ''' Initialise the integrator.

    The integrator requires a manager class implementing the IntegratorManager
    interface. This class executes all the workers in separate threads and
    accumulates the results to expose to the user.

    Params:
      manager The integration manager
      max_procs The number of processors
      mp_method The multiprocessing method

    '''
    self._manager = manager
    self._max_procs = max_procs
    self._mp_method = mp_method

  def integrate(self):
    ''' Do all the integration tasks.

    Returns
      The integration results

    '''
    from time import time
    from libtbx import easy_mp
    start_time = time()
    num_proc = len(self._manager)
    if self._max_procs > 0:
      num_proc = min(num_proc, self._max_procs)
    if num_proc > 1:
      def process_output(result):
        self._manager.accumulate(result[0])
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
        callback=process_output,
        method=self._mp_method,
        preserve_order=True,
        preserve_exception_message=True)
      task_results, output = zip(*task_results)
    else:
      task_results = [task() for task in self._manager.tasks()]
      for result in task_results:
        self._manager.accumulate(result)
    assert(self._manager.finished())
    end_time = time()
    read_time = self._manager.read_time
    extract_time = self._manager.extract_time
    process_time = self._manager.process_time
    total_time = end_time - start_time
    print "Time taken: reading images: %.2f seconds" % read_time
    print "Time taken: extracting pixels: %.2f seconds" % extract_time
    print "Time taken: processing data: %.2f seconds" % process_time
    print "Time taken: total: %.2f seconds" % total_time
    return self._manager.result()


class IntegrationProcessor3D(object):
  ''' A class to do the processing. '''

  def __init__(self, experiments):
    ''' Initialise the processor and set the algorithms. '''
    from dials.framework.registry import Registry
    params = Registry().params()

    # Save the list of experiments
    self._experiments = experiments
    self._n_sigma = params.integration.shoebox.n_sigma
    self._sigma_b = params.integration.shoebox.sigma_b
    self._sigma_m = params.integration.shoebox.sigma_m
    self.time = 0

  def __call__(self, reflections):
    ''' Perform all the processing for integration. '''
    from time import time
    print "Process"
    st = time()

    # Compute the shoebox mask
    reflections.compute_mask(self._experiments,
                             self._n_sigma,
                             self._sigma_b,
                             self._sigma_m)

    # Integrate the reflections
    # reflections.integrate(self._experiments)

    self.time += time() - st
    return reflections

  def compute_mask(self, reflections):
    pass


class IntegrationTask3DExecutorAux(boost.python.injector,
                                   IntegrationTask3DExecutor):
  ''' Add additional methods. '''

  def execute(self, imageset, mask=None):
    ''' Passing in an imageset process all the images. '''
    from dials.model.data import Image
    from time import time
    import sys
    detector = imageset.get_detector()
    frame00 = self.frame0()
    frame10 = self.frame1()
    frame01, frame11 = imageset.get_array_range()
    assert(frame00 == frame01)
    assert(frame10 == frame11)
    self.read_time = 0
    self.extract_time = 0
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
    for i in range(len(imageset)):
      sys.stdout.write(".")
      sys.stdout.flush()
      st = time()
      image = imageset[i]
      self.read_time += time() - st
      if not isinstance(image, tuple):
        image = (image,)
      st = time()
      self.next(Image(image, mask))
      self.extract_time += time() - st
      del image
    sys.stdout.write("\n")
    sys.stdout.flush()
    assert(self.finished())


class IntegrationTask3DMultiExecutor(IntegrationTask3DMultiExecutorBase):
  ''' Add additional methods. '''

  def __init__(self, data, jobs, npanels, callback):
    from dials.array_family import flex

    # Allocate shoeboxes
    data["shoebox"] = flex.shoebox(data["panel"], data["bbox"])
    data["shoebox"].allocate()

    # Set the callback
    self._callback = callback

    # Initialise the base class
    super(IntegrationTask3DMultiExecutor, self).__init__(data, jobs, npanels)

  def execute(self, imageset, mask=None):
    ''' Passing in an imageset process all the images. '''
    from dials.model.data import Image
    from time import time
    import sys
    detector = imageset.get_detector()
    frame00 = self.frame0()
    frame10 = self.frame1()
    frame01, frame11 = imageset.get_array_range()
    assert(frame00 == frame01)
    assert(frame10 == frame11)
    self.read_time = 0
    self.extract_time = 0
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
    for i in range(len(imageset)):
      st = time()
      image = imageset[i]
      self.read_time += time() - st
      if not isinstance(image, tuple):
        image = (image,)
      st = time()
      self.next(Image(image, mask))
      self.extract_time += time() - st
      sys.stdout.write(".")
      sys.stdout.flush()
      del image
    self._callback(self.data())
    del self.data()["shoebox"]
    sys.stdout.write("\n")
    sys.stdout.flush()
    assert(self.finished())


class IntegrationTask3D(IntegrationTask):
  ''' A class to perform a 3D integration task. '''

  def __init__(self, index, experiments, data, jobs):
    ''' Initialise the task. '''
    self._index = index
    self._experiments = experiments
    self._data = data
    self._jobs = jobs
    self._mask = None

  def __call__(self):
    ''' Do the integration. '''
    from dials.array_family import flex
    from scitbx.array_family import shared
    import sys
    print "=" * 80
    print ""
    print "Integrating task %d" % self._index
    process = IntegrationProcessor3D(self._experiments)

    if isinstance(self._jobs, shared.tiny_int_2):
      Executor = IntegrationTask3DExecutor
    else:
      Executor = IntegrationTask3DMultiExecutor
    detectors = self._experiments.detectors()
    assert(len(detectors) == 1)
    npanels = len(detectors[0])
    executor = Executor(self._data, self._jobs, npanels, process)
    imageset = self._experiments[0].imageset
    imageset = imageset[executor.frame0():executor.frame1()]
    executor.execute(imageset)
    result = IntegrationResult(self._index, executor.data())
    result.read_time = executor.read_time
    result.extract_time = executor.extract_time
    result.process_time = process.time
    return result


class IntegrationManager3D(IntegrationManager):
  ''' An class to manage 3D integration. book-keeping '''

  def __init__(self,
               experiments,
               reflections,
               block_size=1,
               max_procs=1):
    ''' Initialise the manager. '''
    from dials.algorithms.integration import IntegrationManager3DExecutor
    from dials.algorithms.integration import IntegrationManager3DMultiExecutor
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
    block_size = block_size / dphi
    if max_procs == 1:
      self._manager = IntegrationManager3DExecutor(
        self._reflections,
        scan.get_array_range(),
        block_size)
    else:
      self._manager = IntegrationManager3DMultiExecutor(
        self._reflections,
        scan.get_array_range(),
        block_size)
    self.read_time = 0
    self.extract_time = 0
    self.process_time = 0

  def task(self, index):
    ''' Get a task. '''
    return IntegrationTask3D(
      index,
      self._experiments,
      self._manager.split(index),
      self._manager.job(index))

  def tasks(self):
    ''' Iterate through the tasks. '''
    for i in range(len(self)):
      yield self.task(i)

  def accumulate(self, result):
    ''' Accumulate the results. '''
    self._manager.accumulate(result.index, result.reflections)
    self.read_time += result.read_time
    self.extract_time += result.extract_time
    self.process_time += result.process_time

  def result(self):
    ''' Return the result. '''
    return self._manager.data()

  def finished(self):
    ''' Return if all tasks have finished. '''
    return self._manager.finished()

  def __len__(self):
    ''' Return the number of tasks. '''
    return len(self._manager)





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

  def __init__(self,
               experiments,
               reflections,
               block_size=1,
               max_procs=1,
               mp_method='multiprocessing'):
    ''' Initialise the manager and the integrator. '''

    # Create the integration manager
    manager = IntegrationManager3D(
      experiments,
      reflections,
      block_size,
      max_procs)

    # Initialise the integrator
    super(Integrator3D, self).__init__(manager, max_procs, mp_method)


class IntegratorFlat2D(Integrator):
  ''' Top level integrator for flat 2D integration. '''

  def __init__(self,
               experiments,
               reflections,
               block_size=1,
               max_procs=1,
               mp_method='multiprocessing'):
    ''' Initialise the manager and the integrator. '''
    raise RuntimeError("Not Implemented")


class Integrator2D(Integrator):
  ''' Top level integrator for 2D integration. '''

  def __init__(self,
               experiments,
               reflections,
               block_size=1,
               max_procs=1,
               mp_method='multiprocessing'):
    ''' Initialise the manager and the integrator. '''
    raise RuntimeError("Not Implemented")


class IntegratorFactory(object):
  ''' A factory for creating integrators. '''

  @staticmethod
  def create(params, experiments, reflections):
    ''' Create the integrator from the input configuration. '''
    from dials.interfaces import IntensityIface
    from dials.interfaces import BackgroundIface
    from dials.interfaces import CentroidIface
    from dials.array_family import flex

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
    IntegratorClass = IntegratorFactory.select_integrator(IntensityAlgorithm)

    # Return an instantiation of the class
    return IntegratorClass(
      experiments=experiments,
      reflections=reflections,
      block_size=params.integration.block.size,
      max_procs=params.integration.mp.max_procs,
      mp_method=params.integration.mp.method)

  @staticmethod
  def select_integrator(cls):
    ''' Select the integrator. '''
    from dials.interfaces import Integration3DMixin
    from dials.interfaces import IntegrationFlat2DMixin
    from dials.interfaces import Integration2DMixin
    if issubclass(cls, Integration3DMixin):
      IntegratorClass = Integrator3D
    elif issubclass(cls, IntegrationFlat2DMixin):
      IntegratorClass = IntegratorFlat2D
    elif issubclass(cls, Integration2DMixin):
      IntegratorClass = Integrator2D
    else:
      raise RuntimeError("Unknown integration type")
    return IntegratorClass
