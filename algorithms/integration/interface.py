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


class IntegrationTask3D(IntegrationTask):

  def __init__(self, index):
    self._index = index

  def __call__(self):
    from time import sleep
    print 'Completing Task %d'% self._index
    sleep(self._index)
    return IntegrationResult(self._index, None)


class IntegrationManager3D(IntegrationManager):

  def __init__(self, experiments, reflections, num_tasks=1, max_overlap=0):
    from dials.array_family import flex
    self._num_tasks = num_tasks
    self._finished = flex.bool(num_tasks, False)

  def task(self, index):
    return IntegrationTask3D(index)

  def tasks(self):
    for i in range(len(self)):
      yield self.task(i)

  def accumulate(self, result):
    assert(result.index < len(self))
    self._finished[result.index] = True

  def result(self):
    assert(self.finished())
    return None

  def finished(self):
    return self._finished.all_eq(True)

  def __len__(self):
    return self._num_tasks

  def _compute_tasks(self, num):
    pass


class Integrator3D(Integrator):
  ''' Top level integrator for 3D integration. '''

  def __init__(self,
               experiments,
               reflections,
               num_tasks=1,
               mp_method='multiprocessing',
               max_procs=None):

    # The 3D Integration manager
    manager = IntegrationManager3D(
      experiments=experiments,
      reflections=reflections,
      num_tasks=num_tasks)

    # Initialise the base integrator
    super(Integrator3D, self).__init__(
      manager=manager,
      mp_method=mp_method,
      max_procs=max_procs)
