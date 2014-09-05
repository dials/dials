from __future__ import division
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
    from libtbx.table_utils import format as table
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
    rows = [
      ["Read time"   , "%.2f seconds" % (self._manager.read_time)   ],
      ["Extract time", "%.2f seconds" % (self._manager.extract_time)],
      ["Process time", "%.2f seconds" % (self._manager.process_time)],
      ["Total time"  , "%.2f seconds" % (end_time - start_time)     ]
    ]
    print table(rows, justify='right', prefix=' ')
    return self._manager.result()


class IntegrationTask3D(IntegrationTask):
  ''' A class to perform a 3D integration task. '''

  def __init__(self, index, experiments, profile_model, data, job, flatten):
    ''' Initialise the task. '''
    self._index = index
    self._experiments = experiments
    self._profile_model = profile_model
    self._data = data
    self._job = job
    self._flatten = flatten
    self._mask = None

  def __call__(self):
    ''' Do the integration. '''
    from dials.array_family import flex
    from time import time
    print "=" * 80
    print ""
    print "Integrating task %d" % self._index

    # Get the sub imageset
    imagesets = self._experiments.imagesets()
    assert(len(imagesets) == 1)
    imageset = imagesets[0]
    frame00, frame01 = self._job
    frame10, frame11 = imageset.get_array_range()
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

    # Create the shoeboxes
    self._data["shoebox"] = flex.shoebox(
      self._data["panel"],
      self._data["bbox"],
      allocate=True)

    # Extract the shoeboxes
    read_time, extract_time = self._data.extract_shoeboxes(imageset, self._mask)

    # Optionally flatten the shoeboxes
    if self._flatten:
      self._data['shoebox'].flatten()

    # Process the data
    st = time()
    self._process()
    process_time = time() - st

    # Delete the shoeboxes
    del self._data["shoebox"]

    # Return the result
    result = IntegrationResult(self._index, self._data)
    result.read_time = read_time
    result.extract_time = extract_time
    result.process_time = process_time
    return result

  def _process(self):
    ''' Process the data. '''
    print "Process"
    self._data.compute_mask(self._experiments, self._profile_model)
    self._data.integrate(self._experiments, self._profile_model)


class IntegrationManager3D(IntegrationManager):
  ''' An class to manage 3D integration. book-keeping '''

  def __init__(self,
               experiments,
               profile_model,
               reflections,
               block_size=1,
               max_procs=1,
               flatten=False):
    ''' Initialise the manager. '''
    from dials.algorithms.integration import IntegrationManager3DExecutor
    from dials.array_family import flex
    imagesets = experiments.imagesets()
    detectors = experiments.detectors()
    scans = experiments.scans()
    assert(len(experiments) == len(profile_model))
    assert(len(imagesets) == 1)
    assert(len(scans) == 1)
    assert(len(detectors) == 1)
    imageset = imagesets[0]
    scan = scans[0]
    detector = detectors[0]
    assert(len(imageset) == len(scan))
    self._flatten = flatten
    self._experiments = experiments
    self._profile_model = profile_model
    self._reflections = reflections
    phi0, dphi = scan.get_oscillation()
    block_size = block_size / dphi
    self._manager = IntegrationManager3DExecutor(
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
      self._profile_model,
      self._manager.split(index),
      self._manager.job(index),
      self._flatten)

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
    if self.finished():
      self._post_process(self._manager.data())

  def result(self):
    ''' Return the result. '''
    return self._manager.data()

  def finished(self):
    ''' Return if all tasks have finished. '''
    return self._manager.finished()

  def __len__(self):
    ''' Return the number of tasks. '''
    return len(self._manager)

  def _post_process(self, data):
    ''' Do some post processing. '''
    data.compute_corrections(self._experiments)


class Integrator3D(Integrator):
  ''' Top level integrator for 3D integration. '''

  def __init__(self,
               experiments,
               profile_model,
               reflections,
               block_size=1,
               max_procs=1,
               mp_method='multiprocessing'):
    ''' Initialise the manager and the integrator. '''

    # Create the integration manager
    manager = IntegrationManager3D(
      experiments,
      profile_model,
      reflections,
      block_size,
      max_procs)

    # Initialise the integrator
    super(Integrator3D, self).__init__(manager, max_procs, mp_method)


class IntegratorFlat2D(Integrator):
  ''' Top level integrator for flat 2D integration. '''

  def __init__(self,
               experiments,
               profile_model,
               reflections,
               block_size=1,
               max_procs=1,
               mp_method='multiprocessing'):
    ''' Initialise the manager and the integrator. '''

    # Create the integration manager
    manager = IntegrationManager3D(
      experiments,
      profile_model,
      reflections,
      block_size,
      max_procs,
      flatten=True)

    # Initialise the integrator
    super(Integrator3D, self).__init__(manager, max_procs, mp_method)


class Integrator2D(Integrator):
  ''' Top level integrator for 2D integration. '''

  def __init__(self,
               experiments,
               profile_model,
               reflections,
               block_size=1,
               max_procs=1,
               mp_method='multiprocessing'):
    ''' Initialise the manager and the integrator. '''
    raise RuntimeError("Not Implemented")


class IntegratorFactory(object):
  ''' A factory for creating integrators. '''

  @staticmethod
  def create(params, experiments, profile_model, reflections):
    ''' Create the integrator from the input configuration. '''
    from dials.interfaces import IntensityIface
    from dials.interfaces import BackgroundIface
    from dials.interfaces import CentroidIface
    from dials.array_family import flex

    # Check the input
    assert(len(experiments) == len(profile_model))

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
      profile_model=profile_model,
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
