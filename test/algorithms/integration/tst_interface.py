
from __future__ import division

class TestIntegrationManagerExecutor(object):

  def __init__(self):
    from dials.array_family import flex

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
    from dials.algorithms.integration import IntegrationManagerExecutor
    from dials.array_family import flex

    # Create the executor
    executor = IntegrationManagerExecutor(
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
    from dials.algorithms.profile_model.profile_model import ProfileModelList
    from dials.algorithms.profile_model.profile_model import ProfileModel
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
    profile_model = ProfileModelList()
    profile_model.append(ProfileModel(
      n_sigma=3,
      sigma_b=0.024*pi/180.0,
      sigma_m=0.044*pi/180.0))
    self.nproc = nproc

    rlist = flex.reflection_table.from_predictions(exlist[0])
    rlist['id'] = flex.size_t(len(rlist), 0)
    rlist.compute_bbox(exlist, profile_model)
    rlist.compute_zeta_multi(exlist)
    rlist.compute_d(exlist)
    self.rlist = rlist
    self.exlist = exlist
    self.profile_model = profile_model

  def run(self):
    from dials.algorithms.integration.interface import Integrator3D

    integrator = Integrator3D(
      self.exlist,
      self.profile_model,
      self.rlist,
      block_size=5,
      nproc=self.nproc)
    result = integrator.integrate()

    print 'OK'


class TestSummation(object):

  def __init__(self):
    from dxtbx.model.experiment.experiment_list import ExperimentListFactory
    from dials.algorithms.profile_model.profile_model import ProfileModelList
    from dials.algorithms.profile_model.profile_model import ProfileModel
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
    profile_model = ProfileModelList()
    profile_model.append(ProfileModel(
      n_sigma=3,
      sigma_b=0.024*pi/180.0,
      sigma_m=0.044*pi/180.0))

    rlist = flex.reflection_table.from_predictions(exlist[0])
    rlist['id'] = flex.size_t(len(rlist), 0)
    self.rlist = rlist
    self.exlist = exlist
    self.profile_model = profile_model

  def run(self):
    from libtbx.test_utils import approx_equal
    from dials.array_family import flex

    def approx_equal_dict(a, b, k):
      return approx_equal(a[k], b[k])

    # Do summation by all different methods
    result1 = self.integrate("3d")
    # result2 = self.integrate("flat3d")
    result3 = self.integrate("2d")
    result4 = self.integrate("single2d")
    assert(len(result1) >= len(self.rlist))
    # assert(len(result2) >= len(self.rlist))
    assert(len(result3) >= len(self.rlist))
    assert(len(result4) >= len(self.rlist))

    # Ensure we get equivalent results

    # result3 and result4 should be the same
    assert(len(result3) == len(result4))
    for r3, r4 in zip(result3, result4):
      assert(r3['partial_id'] == r4['partial_id'])
      assert(r3['bbox'] == r4['bbox'])
      assert(r3['entering'] == r4['entering'])
      assert(r3['flags'] == r4['flags'])
      assert(r3['id'] == r4['id'])
      assert(r3['miller_index'] == r4['miller_index'])
      assert(r3['panel'] == r4['panel'])
      assert(approx_equal_dict(r3, r4, 'd'))
      assert(approx_equal_dict(r3, r4, 'intensity.sum.value'))
      assert(approx_equal_dict(r3, r4, 'intensity.sum.variance'))
      assert(approx_equal_dict(r3, r4, 'lp'))
      assert(approx_equal_dict(r3, r4, 'partiality'))
      assert(approx_equal_dict(r3, r4, 's1'))
      assert(approx_equal_dict(r3, r4, 'xyzcal.mm'))
      assert(approx_equal_dict(r3, r4, 'xyzcal.px'))
      assert(approx_equal_dict(r3, r4, 'xyzobs.px.value'))
      assert(approx_equal_dict(r3, r4, 'xyzobs.px.variance'))
      assert(approx_equal_dict(r3, r4, 'zeta'))

    # result3 should add up to result1
    assert(len(result3) >= len(result1))
    expected1 = self.rlist.copy()
    expected1['intensity.sum.value'] = flex.double(len(self.rlist), 0)
    expected1['intensity.sum.variance'] = flex.double(len(self.rlist), 0)
    for r1 in result1:
      pid = r1['partial_id']
      r2 = expected1[pid]
      assert(r1['entering'] == r2['entering'])
      assert(r1['id'] == r2['id'])
      assert(r1['miller_index'] == r2['miller_index'])
      assert(r1['panel'] == r2['panel'])
      assert(approx_equal_dict(r1, r2, 's1'))
      assert(approx_equal_dict(r1, r2, 'xyzcal.mm'))
      assert(approx_equal_dict(r1, r2, 'xyzcal.px'))
      expected1['intensity.sum.value'][pid] += r1['intensity.sum.value']
      expected1['intensity.sum.variance'][pid] += r1['intensity.sum.variance']
    expected3 = self.rlist.copy()
    expected3['intensity.sum.value'] = flex.double(len(self.rlist), 0)
    expected3['intensity.sum.variance'] = flex.double(len(self.rlist), 0)
    for r1 in result3:
      pid = r1['partial_id']
      r2 = expected3[pid]
      assert(r1['entering'] == r2['entering'])
      assert(r1['id'] == r2['id'])
      assert(r1['miller_index'] == r2['miller_index'])
      assert(r1['panel'] == r2['panel'])
      assert(approx_equal_dict(r1, r2, 's1'))
      assert(approx_equal_dict(r1, r2, 'xyzcal.mm'))
      assert(approx_equal_dict(r1, r2, 'xyzcal.px'))
      expected3['intensity.sum.value'][pid] += r1['intensity.sum.value']
      expected3['intensity.sum.variance'][pid] += r1['intensity.sum.variance']
    for r1, r3, in zip(expected1, expected3):
      assert(approx_equal_dict(r1, r3, 'intensity.sum.value'))
      assert(approx_equal_dict(r1, r3, 'intensity.sum.variance'))


    print 'OK'


  def integrate(self, integrator_type):
    from dials.algorithms.integration.interface import IntegratorFactory
    from dials.algorithms.integration.interface import phil_scope as master_phil_scope
    from libtbx.phil import parse
    import sys
    import StringIO

    rlist = self.rlist.copy()

    output = StringIO.StringIO()
    stdout = sys.stdout
    sys.stdout = output

    try:
      phil_scope = parse('''
        integration.background.algorithm=null
        integration.intensity.algorithm=sum
        integration.intensity.sum.integrator=%s
        integration.block.size=0.5
      ''' % integrator_type)

      params = master_phil_scope.fetch(source=phil_scope).extract()

      integrator = IntegratorFactory.create(
        params,
        self.exlist,
        self.profile_model,
        rlist)

      result = integrator.integrate()
    except Exception:
      print output
      raise

    sys.stdout = stdout

    return result


class Test(object):

  def __init__(self):
    self.test1 = TestIntegrationManagerExecutor()
    self.test2 = TestIntegrator3D(nproc=1)
    self.test3 = TestIntegrator3D(nproc=2)
    self.test4 = TestSummation()

  def run(self):
    # self.test1.run()
    # self.test2.run()
    # self.test3.run()
    self.test4.run()

if __name__ == '__main__':
  test = Test()
  test.run()
