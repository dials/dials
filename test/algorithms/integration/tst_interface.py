
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
    self.nrefl = 10000
    # self.nrefl = 100000

    self.jobs = shared.tiny_int_2([
      (0, 16),
      (8, 24),
      (16, 32),
      (24, 40)])

    self.jobs = shared.tiny_int_2([
      (0, 40),
      (20, 60),
      (40, 80),
      (80, 100)])

    for i, j in enumerate(self.jobs):
      self.append_reflections(j, i)

    ind1 = flex.size_t(range(self.nrefl))
    ind2 = ind1 + self.nrefl
    ind3 = ind2 + self.nrefl
    ind4 = ind3 + self.nrefl
    self.indices = ind1
    self.indices.extend(ind2)
    self.indices.extend(ind3)
    self.indices.extend(ind4)
    off1 = self.nrefl
    off2 = off1 + self.nrefl
    off3 = off2 + self.nrefl
    off4 = off3 + self.nrefl
    self.offset = flex.size_t([0, off1, off2, off3, off4])

    # indices = list(range(len(self.reflections)))
    # shuffle(indices)
    # self.reflections.reorder(flex.size_t(indices))

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
    from dials.algorithms.integration import IntegrationTask3DSpec
    from dials.array_family import flex
    from dials.model.data import Image
    from time import time

    # The processing callback
    class Callback(object):
      def __init__(self, nrefl):
        self.ncallback = 0
        self.nrefl = nrefl
      def __call__(self, reflections):
        print reflections
        assert(len(reflections) == self.nrefl)
        for sbox in reflections['shoebox']:
          bbox = sbox.bbox
          v1 = bbox[4]+1
          for z in range(bbox[5]-bbox[4]):
            assert(sbox.data[z:z+1,:,:].all_eq(v1+z))
        self.ncallback += 1
        return reflections
    callback = Callback(self.nrefl)

    print len(self.reflections)

    # Create the task specification
    st = time()
    spec = IntegrationTask3DSpec(
      self.reflections,
      self.npanels,
      self.jobs,
      self.offset,
      self.indices)

    # Initialise the executor
    executor = IntegrationTask3DExecutor(spec, callback)
    print time() - st

    # Check the initial state is correct
    assert(executor.frame0() == 0)
    assert(executor.frame1() == 100)
    assert(executor.nframes() == 100)

    # The data and mask
    data = flex.int(flex.grid(self.height,self.width), 1)
    mask = flex.bool(flex.grid(self.height,self.width), True)

    # Loop through images
    st = time()
    for i in range(100):
      data += 1
      executor.next(Image((data, data), (mask, mask)))
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
