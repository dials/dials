
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

    # self.jobs = shared.tiny_int_2([
    #   (0, 16),
    #   (8, 24),
    #   (16, 32),
    #   (24, 40)])

    self.jobs = shared.tiny_int_2([
      (0, 40),
      (20, 60),
      (40, 80),
      (60, 100)])

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
      zs = randint(1, (zrange[1] - zrange[0]) // 2)
      x1 = x0 + randint(1, 10)
      y1 = y0 + randint(1, 10)
      z0 = (zrange[1] + zrange[0]) // 2 - zs // 2
      z1 = z0 + zs
      assert(z1 > z0)
      assert(z0 >= zrange[0] and z1 <= zrange[1])
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
    from dials.array_family import flex
    from dials.model.data import Image
    from time import time

    # The processing callback
    class Callback(object):
      def __init__(self, nrefl):
        self.ncallback = 0
        self.nrefl = nrefl
      def __call__(self, reflections):
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

    # Initialise the executor
    st = time()
    executor = IntegrationTask3DExecutor(
      self.reflections,
      self.jobs,
      self.npanels,
      callback)
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

class TestIntegrationTask3DExecutorMulti(object):

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

    # self.jobs = shared.tiny_int_2([
    #   (0, 16),
    #   (8, 24),
    #   (16, 32),
    #   (24, 40)])

    self.jobs = shared.tiny_int_2([
      (0, 40),
      (20, 60),
      (40, 80),
      (60, 100)])

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
      zs = randint(1, (zrange[1] - zrange[0]) // 2)
      x1 = x0 + randint(1, 10)
      y1 = y0 + randint(1, 10)
      z0 = (zrange[1] + zrange[0]) // 2 - zs // 2
      z1 = z0 + zs
      assert(x1 > x0)
      assert(y1 > y0)
      assert(z1 > z0)
      assert(z0 >= zrange[0] and z1 <= zrange[1])
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

    from dials.algorithms.integration import IntegrationTask3DMultiExecutorBase
    from dials.array_family import flex
    from dials.model.data import Image
    from time import time

    # The processing callback
    class Callback(object):
      def __init__(self, nrefl):
        self.ncallback = 0
        self.nrefl = nrefl
      def __call__(self, reflections):
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

    # Initialise the executor
    for j, job in enumerate(self.jobs):
      indices = flex.size_t(range(10000)) + 10000*j
      reflections = self.reflections.select(indices)
      reflections["shoebox"] = flex.shoebox(
        reflections["panel"],
        reflections["bbox"])
      reflections["shoebox"].allocate()
      executor = IntegrationTask3DMultiExecutorBase(
        reflections,
        job,
        self.npanels)

      # Check the initial state is correct
      assert(executor.frame0() == job[0])
      assert(executor.frame1() == job[1])
      assert(executor.nframes() == job[1] - job[0])

      # The data and mask
      data = flex.int(flex.grid(self.height,self.width), job[0]+1)
      mask = flex.bool(flex.grid(self.height,self.width), True)

      # Loop through images
      st = time()
      for i in range(job[0], job[1]):
        executor.next(Image((data, data), (mask, mask)))
        data += 1
      # print time() - st
      assert(executor.finished())
      callback(executor.data())

      del reflections["shoebox"]

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
    self.nrefl = 10000
    self.array_range = (0, 130)
    self.block_size = 20

    from random import randint, seed, choice
    seed(0)
    self.expected = []
    self.processed = []
    for i in range(self.nrefl):
      x0 = randint(0, self.width-10)
      y0 = randint(0, self.height-10)
      zs = randint(2, 9)
      x1 = x0 + randint(2, 10)
      y1 = y0 + randint(2, 10)
      for k, j in enumerate([10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 110, 120]):
        m = k + i * 13
        pos = choice(["left", "right", "centre"])
        if pos == 'left':
          z0 = j - zs
          z1 = j
        elif pos == 'right':
          z0 = j
          z1 = j + zs
        else:
          z0 = j - zs // 2
          z1 = j + zs // 2
        bbox = (x0, x1, y0, y1, z0, z1)
        self.reflections.append({
          "panel" : randint(0,1),
          "bbox" : bbox,
          "flags" : flex.reflection_table.flags.reference_spot
        })
        self.expected.append(m)
        self.processed.append(m)

      # Add reflection to ignore
      zc = choice([10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 110, 120])
      z0 = zc - 11
      z1 = zc + 11
      bbox = (x0, x1, y0, y1, z0, z1)
      self.reflections.append({
        "panel" : randint(0,1),
        "bbox" : bbox,
        "flags" : flex.reflection_table.flags.reference_spot
      })


  def run(self):
    from dials.algorithms.integration import IntegrationManager3DExecutor
    from dials.array_family import flex

    # Create the executor
    executor = IntegrationManager3DExecutor(
      self.reflections,
      self.array_range,
      self.block_size)

    # Ensure the tasks make sense
    jobs = executor.job(0)
    assert(len(executor) == 1)
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

    # Get the task specs
    data = executor.split(0)
    ignored = executor.ignored()
    assert(len(data) == len(self.expected))
    assert(len(executor.ignored()) == self.nrefl)
    assert(len(data) + len(executor.ignored()) == len(self.reflections))

    # Add some results
    data["data"] = flex.double(len(data), 1)

    # Accumulate the data again
    assert(not executor.finished())
    executor.accumulate(0, data)
    assert(executor.finished())

    # Get results and check they're as expected
    data = executor.data()
    result = data["data"]
    bbox = data["bbox"]
    for i in range(len(self.processed)):
      assert(result[self.processed[i]] == 1)

    # Test passed
    print 'OK'


class TestIntegrationManager3DExecutorMulti(object):

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
    self.array_range = (0, 130)
    self.block_size = 20

    from random import randint, seed, choice
    seed(0)
    self.expected = [[] for i in range(12)]
    self.processed = [[] for i in range(12)]
    for i in range(self.nrefl):
      x0 = randint(0, self.width-10)
      y0 = randint(0, self.height-10)
      zs = randint(2, 9)
      x1 = x0 + randint(2, 10)
      y1 = y0 + randint(2, 10)
      for k, j in enumerate([10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 110, 120]):
        m = k + i * 13
        pos = choice(["left", "right", "centre"])
        if pos == 'left':
          z0 = j - zs
          z1 = j
          if k > 0:
            self.expected[k-1].append(m)
        elif pos == 'right':
          z0 = j
          z1 = j + zs
          if k < 11:
            self.expected[k+1].append(m)
        else:
          z0 = j - zs // 2
          z1 = j + zs // 2
        bbox = (x0, x1, y0, y1, z0, z1)
        self.reflections.append({
          "panel" : randint(0,1),
          "bbox" : bbox,
          "flags" : flex.reflection_table.flags.reference_spot
        })
        self.expected[k].append(m)
        self.processed[k].append(m)

      # Add reflection to ignore
      zc = choice([10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 110, 120])
      z0 = zc - 11
      z1 = zc + 11
      bbox = (x0, x1, y0, y1, z0, z1)
      self.reflections.append({
        "panel" : randint(0,1),
        "bbox" : bbox,
        "flags" : flex.reflection_table.flags.reference_spot
      })


  def run(self):
    from dials.algorithms.integration import IntegrationManager3DMultiExecutor
    from dials.array_family import flex

    # Create the executor
    executor = IntegrationManager3DMultiExecutor(
      self.reflections,
      self.array_range,
      self.block_size)

    # Ensure the tasks make sense
    jobs = [executor.job(i) for i in range(len(executor))]
    assert(len(executor) == 12)
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

    # Get the task specs
    data0 = executor.split(0)
    data1 = executor.split(1)
    data2 = executor.split(2)
    data3 = executor.split(3)
    data4 = executor.split(4)
    data5 = executor.split(5)
    data6 = executor.split(6)
    data7 = executor.split(7)
    data8 = executor.split(8)
    data9 = executor.split(9)
    data10 = executor.split(10)
    data11 = executor.split(11)
    ignored = executor.ignored()
    assert(len(data0) == len(self.expected[0]))
    assert(len(data1) == len(self.expected[1]))
    assert(len(data2) == len(self.expected[2]))
    assert(len(data3) == len(self.expected[3]))
    assert(len(data4) == len(self.expected[4]))
    assert(len(data5) == len(self.expected[5]))
    assert(len(data6) == len(self.expected[6]))
    assert(len(data7) == len(self.expected[7]))
    assert(len(data8) == len(self.expected[8]))
    assert(len(data9) == len(self.expected[9]))
    assert(len(data10) == len(self.expected[10]))
    assert(len(data11) == len(self.expected[11]))
    assert(len(executor.ignored()) == self.nrefl)

    # Add some results
    data0["data"] = flex.double(len(data0), 1)
    data1["data"] = flex.double(len(data1), 2)
    data2["data"] = flex.double(len(data2), 3)
    data3["data"] = flex.double(len(data3), 4)
    data4["data"] = flex.double(len(data4), 5)
    data5["data"] = flex.double(len(data5), 6)
    data6["data"] = flex.double(len(data6), 7)
    data7["data"] = flex.double(len(data7), 8)
    data8["data"] = flex.double(len(data8), 9)
    data9["data"] = flex.double(len(data9), 10)
    data10["data"] = flex.double(len(data10), 11)
    data11["data"] = flex.double(len(data11), 12)

    # Accumulate the data again
    assert(not executor.finished())
    executor.accumulate(0, data0)
    executor.accumulate(1, data1)
    executor.accumulate(2, data2)
    executor.accumulate(3, data3)
    executor.accumulate(4, data4)
    executor.accumulate(5, data5)
    executor.accumulate(6, data6)
    executor.accumulate(7, data7)
    executor.accumulate(8, data8)
    executor.accumulate(9, data9)
    executor.accumulate(10, data10)
    executor.accumulate(11, data11)
    assert(executor.finished())

    # Get results and check they're as expected
    data = executor.data()
    result = data["data"]
    bbox = data["bbox"]
    for i in range(len(self.processed)):
      for j in range(len(self.processed[i])):
        assert(result[self.processed[i][j]] == i+1)

    # Test passed
    print 'OK'


class TestIntegrator3D(object):

  def __init__(self, nproc):
    from dxtbx.model.experiment.experiment_list import ExperimentListFactory
    import libtbx.load_env
    from dials.array_family import flex
    from os.path import join
    from math import pi
    try:
      dials_regression = libtbx.env.dist_path('dials_regression')
    except KeyError, e:
      print 'FAIL: dials_regression not configured'
      exit(0)


    path = join(dials_regression, "centroid_test_data", "experiments.json")

    exlist = ExperimentListFactory.from_json_file(path)

    self.nproc = nproc

    rlist = flex.reflection_table.from_predictions(exlist[0])
    rlist['id'] = flex.size_t(len(rlist), 0)
    rlist.compute_bbox(exlist[0], nsigma=3, sigma_d=0.024*pi/180,
                       sigma_m=0.044*pi/180)
    rlist.compute_zeta_multi(exlist)
    rlist.compute_d(exlist)
    self.rlist = rlist
    self.exlist = exlist

  def run(self):
    from dials.algorithms.integration.interface import Integrator3D
    from dials.algorithms.integration.interface import phil_scope
    from libtbx import phil

    user_phil = phil.parse('''
      mp.max_procs = %d
      block.size=5
      filter.ice_rings.filter=False
    ''' % self.nproc)

    params = phil_scope.fetch(source=user_phil).extract()

    integrator = Integrator3D(self.exlist, self.rlist, params)
    result = integrator.integrate()

    print 'OK'


class Test(object):

  def __init__(self):
    self.test1 = TestIntegrationTask3DExecutor()
    self.test2 = TestIntegrationTask3DExecutorMulti()
    self.test3 = TestIntegrationManager3DExecutor()
    self.test4 = TestIntegrationManager3DExecutorMulti()
    self.test5 = TestIntegrator3D(nproc=1)
    self.test6 = TestIntegrator3D(nproc=2)

  def run(self):
    self.test1.run()
    self.test2.run()
    self.test3.run()
    self.test4.run()
    self.test5.run()
    self.test6.run()

if __name__ == '__main__':
  test = Test()
  test.run()
