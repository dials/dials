
from __future__ import division

class TestIntegrationTask3DExecutor(object):

  def __init__(self):
    from dials.array_family import flex
    from scitbx.array_family import shared
    from random import shuffle

    self.reflections = flex.reflection_table()
    self.reflections['panel'] = flex.size_t()
    self.reflections['bbox'] = flex.int6()
    self.reflections['miller_index'] = flex.miller_index()
    self.reflections['s1'] = flex.vec3_double()
    self.reflections['xyzcal.px'] = flex.vec3_double()
    self.reflections['xyzcal.mm'] = flex.vec3_double()
    self.reflections['entering'] = flex.bool()
    self.reflections['id'] = flex.size_t()
    self.reflections["flags"] = flex.size_t()

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
      self.append_reflections(j)

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
    self.mask = flex.bool(len(self.indices), True)

    # indices = list(range(len(self.reflections)))
    # shuffle(indices)
    # self.reflections.reorder(flex.size_t(indices))

  def append_reflections(self, zrange):
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
        # print reflections, len(reflections)
        assert(len(reflections) == self.nrefl)
        for sbox in reflections['shoebox']:
          assert(sbox.is_consistent())
          bbox = sbox.bbox
          v1 = bbox[4]+1
          for z in range(bbox[5]-bbox[4]):
            assert(sbox.data[z:z+1,:,:].all_eq(v1+z))
        self.ncallback += 1
        return reflections
    callback = Callback(self.nrefl)

    # print len(self.reflections)

    # Create the task specification
    st = time()
    spec = IntegrationTask3DSpec(
      self.reflections,
      self.npanels,
      self.jobs,
      self.offset,
      self.indices,
      self.mask)
    # print time() - st

    # Initialise the executor
    st = time()
    executor = IntegrationTask3DExecutor(spec, callback)
    # print time() - st

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
      executor.next(Image((data, data), (mask, mask)))
      data += 1
    # print time() - st
    assert(executor.finished())
    assert(callback.ncallback == 4)
    print 'OK'


class TestIntegrationManager3DExecutor(object):

  def __init__(self):
    from dials.array_family import flex
    from scitbx.array_family import shared
    from random import shuffle

    self.reflections = flex.reflection_table()
    self.reflections['panel'] = flex.size_t()
    self.reflections['bbox'] = flex.int6()
    self.reflections['miller_index'] = flex.miller_index()
    self.reflections['s1'] = flex.vec3_double()
    self.reflections['xyzcal.px'] = flex.vec3_double()
    self.reflections['xyzcal.mm'] = flex.vec3_double()
    self.reflections['entering'] = flex.bool()
    self.reflections['id'] = flex.size_t()
    self.reflections["flags"] = flex.size_t()

    self.npanels = 2
    self.width = 1000
    self.height = 1000
    self.nrefl = 100000
    self.array_range = (0, 130)
    self.block_size = 20
    self.num_tasks = 4

    from random import randint, seed
    seed(0)
    for i in range(self.nrefl):
      x0 = randint(0, self.width-10)
      y0 = randint(0, self.height-10)
      z0 = randint(self.array_range[0], self.array_range[1]-1)
      x1 = x0 + randint(1, 10)
      y1 = y0 + randint(1, 10)
      z1 = randint(z0+1, min(self.array_range[1], z0+10))
      bbox = (x0, x1, y0, y1, z0, z1)
      self.reflections.append({
        "panel" : randint(0,1),
        "bbox" : bbox,
      })

  def run(self):
    from dials.algorithms.integration import IntegrationManager3DExecutor

    # Create the executor
    executor = IntegrationManager3DExecutor(
      self.reflections,
      self.array_range,
      self.block_size,
      self.num_tasks,
      self.npanels)

    # Ensure the tasks make sense
    jobs = executor.jobs()
    assert(len(executor) == 4)
    assert(not executor.finished())
    assert(len(jobs) == 12)
    assert(jobs[0] == (0, 20))
    assert(jobs[1] == (10, 30))
    assert(jobs[2] == (20, 40))
    assert(jobs[3] == (30, 50))
    assert(jobs[4] == (40, 60))
    assert(jobs[5] == (50, 70))
    assert(jobs[6] == (60, 80))
    assert(jobs[7] == (70, 90))
    assert(jobs[8] == (80, 100))
    assert(jobs[9] == (90, 110))
    assert(jobs[10] == (100, 120))
    assert(jobs[11] == (110, 130))
    assert(executor.task(0) == (0, 3))
    assert(executor.task(1) == (3, 6))
    assert(executor.task(2) == (6, 9))
    assert(executor.task(3) == (9, 12))

    # Get the task specs
    spec0 = executor.split(0)
    spec1 = executor.split(1)
    spec2 = executor.split(2)
    spec3 = executor.split(3)
    assert(spec0.npanels() == 2)
    assert(spec1.npanels() == 2)
    assert(spec2.npanels() == 2)
    assert(spec3.npanels() == 2)
    assert(spec0.njobs() == 3)
    assert(spec1.njobs() == 3)
    assert(spec2.njobs() == 3)
    assert(spec3.njobs() == 3)
    assert(spec0.frame0() == 0)
    assert(spec0.frame1() == 40)
    assert(spec0.nframes() == 40)
    assert(spec1.frame0() == 30)
    assert(spec1.frame1() == 70)
    assert(spec1.nframes() == 40)
    assert(spec2.frame0() == 60)
    assert(spec2.frame1() == 100)
    assert(spec2.nframes() == 40)
    assert(spec3.frame0() == 90)
    assert(spec3.frame1() == 130)
    assert(spec3.nframes() == 40)
    assert(spec0.job(0) == (0, 20))
    assert(spec0.job(1) == (10, 30))
    assert(spec0.job(2) == (20, 40))
    assert(spec1.job(0) == (30, 50))
    assert(spec1.job(1) == (40, 60))
    assert(spec1.job(2) == (50, 70))
    assert(spec2.job(0) == (60, 80))
    assert(spec2.job(1) == (70, 90))
    assert(spec2.job(2) == (80, 100))
    assert(spec3.job(0) == (90, 110))
    assert(spec3.job(1) == (100, 120))
    assert(spec3.job(2) == (110, 130))

    # Accumulate the data again
    executor.accumulate(0, spec0.data())
    executor.accumulate(1, spec1.data())
    executor.accumulate(2, spec2.data())
    executor.accumulate(3, spec3.data())

    # Test passed
    print 'OK'


class Test(object):

  def __init__(self):
    # self.test1 = TestIntegrationTask3DExecutor()
    self.test2 = TestIntegrationManager3DExecutor()

  def run(self):
    # self.test1.run()
    self.test2.run()

if __name__ == '__main__':
  test = Test()
  test.run()
