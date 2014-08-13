from __future__ import division
import abc

class IntegrationResult(object):
  
  def __init__(self, index):
    self._index = index


class IntegrationTask(object):
  
  __metaclass__ = abc.ABCMeta
  
  @abc.abstractmethod
  def complete(self):
    pass

  @abc.abstractmethod
  def result(self):
    pass


class IntegrationManager(object):
  
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
  def __len__(self):
    pass


class IntegrationTask3D(object):

  def __init__(self, index):
    self._index = index

  def complete(self):
    from time import sleep
    print 'Completing Task ', self._index
    sleep(2)

  def result(self):
    return IntegrationResult(self._index)


class IntegrationManager3D(object):

  def __init__(self, experiments, reflections, ntasks=1):
    self._ntasks = ntasks

  def task(self, index):
    return IntegrationTask3D(index)

  def tasks(self):
    for i in range(len(self)):
      yield self.task(i)

  def accumulate(self, result):
    pass

  def result(self):
    return None

  def __len__(self):
    return self._ntasks


class Integrator(object):

  def __init__(self, manager, mp_method='multiprocessing', max_proc=None):
    self._mp_method = mp_method
    self._manager = manager
    self._max_proc = max_proc

  def integrate(self):
    from libtbx import easy_mp
    num_proc = len(self._manager)
    if self._max_proc is not None:
      num_proc = min(num_proc, self._max_proc)
    if num_proc > 1:
      task_results = easy_mp.parallel_map(
        func=lambda task: task.complete(),
        iterable=list(self._manager.tasks()),
        processes=num_proc,
        method=self._mp_method,
        preserve_order=True)
    else:
      task_results = [self._manager.task(0).complete()]
    for result in task_results:
      self._manager.accumulate(result)
    return self._manager.result()
