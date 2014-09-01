
from __future__ import division

class TestIntegrationTask3DExecutor(object):

  def __init__(self):
    from dials.array_family import flex
    from scitbx.array_family import shared
    from random import shuffle

    self.reflections = flex.reflection_table()
    self.reflections['panel'] = flex.size_t()
    self.reflections['bbox'] = flex.int6()
    self.reflections['job_id'] = flex.size_t()

    self.npanels = 2
    self.width = 1000
    self.height = 1000
    self.nrefl = 20000

    self.jobs = shared.tiny_int_2([
      (0, 4),
      (2, 6),
      (4, 8),
      (6, 10)])

    for i, j in enumerate(self.jobs):
      self.append_reflections(j, i)

    indices = list(range(len(self.reflections)))
    shuffle(indices)
    self.reflections.reorder(flex.size_t(indices))

    self.images = [self.create_image(i) for i in range(10)]

  def append_reflections(self, zrange, job_id):
    from random import randint, seed
    seed(0)
    for i in range(self.nrefl):
      x0 = randint(0, self.width-10)
      y0 = randint(0, self.height-10)
      z0 = randint(zrange[0], zrange[1]-1)
      x1 = x0 + randint(1, 10)
      y1 = y0 + randint(1, 10)
      z1 = randint(z0+1, zrange[1])
      bbox = (x0, x1, y0, y1, z0, z1)
      self.reflections.append({
        "panel" : randint(0,1),
        "bbox" : bbox,
        "job_id" : job_id
      })

  def create_image(self, value):
    from dials.array_family import flex
    from dials.model.data import Image
    data = flex.int(flex.grid(self.height,self.width), value)
    mask = flex.bool(flex.grid(self.height,self.width), True)
    return Image((data, data), (mask, mask))

  def run(self):

    from dials.algorithms.integration import IntegrationTask3DExecutor
    from time import time

    # The processing callback
    class Callback(object):
      def __init__(self, nrefl):
        self.ncallback = 0
        self.nrefl = nrefl
      def __call__(self, reflections):
        print reflections
        assert(len(reflections) == self.nrefl)
        assert(reflections['job_id'].all_eq(self.ncallback))
        self.ncallback += 1
        return reflections
    callback = Callback(self.nrefl)

    # Initialise the executor
    st = time()
    executor = IntegrationTask3DExecutor(
      self.reflections,
      self.jobs,
      self.npanels,
      callback)
    print time() - st

    # Check the initial state is correct
    assert(executor.frame0() == 0)
    assert(executor.frame1() == 10)
    assert(executor.nframes() == 10)

    # Loop through images
    st = time()
    for image in self.images:
      executor.next(image)
    assert(executor.finished())
    assert(callback.ncallback == 4)
    print time() - st


class Test(object):

  def __init__(self):
    self.test1 = TestIntegrationTask3DExecutor()

  def run(self):
    self.test1.run()

if __name__ == '__main__':
  test = Test()
  test.run()
