from __future__ import absolute_import, division, print_function

class Test(object):

  def __init__(self):
    from dials.model.data import Observation
    self.obs = Observation()
    self.obs.panel = 1
    self.obs.centroid.px.position = (100, 200, 300)
    self.obs.centroid.px.variance = (10, 20, 30)
    self.obs.centroid.px.std_err_sq = (1, 2, 3)
    self.obs.centroid.mm.position = (100.1, 200.1, 300.1)
    self.obs.centroid.mm.variance = (10.1, 20.1, 30.1)
    self.obs.centroid.mm.std_err_sq = (1.1, 2.1, 3.1)
    self.obs.intensity.observed.value = 1234
    self.obs.intensity.observed.variance = 42
    self.obs.intensity.corrected.value = 1234.1
    self.obs.intensity.corrected.variance = 42.1

  def run(self):
    self.tst_data()
    self.tst_equality()

  def tst_data(self):
    eps = 1e-7

    assert(self.obs.panel == 1)

    assert(abs(self.obs.centroid.px.position[0] - 100) < eps)
    assert(abs(self.obs.centroid.px.position[1] - 200) < eps)
    assert(abs(self.obs.centroid.px.position[2] - 300) < eps)
    assert(abs(self.obs.centroid.px.variance[0] - 10) < eps)
    assert(abs(self.obs.centroid.px.variance[1] - 20) < eps)
    assert(abs(self.obs.centroid.px.variance[2] - 30) < eps)
    assert(abs(self.obs.centroid.px.std_err_sq[0] - 1) < eps)
    assert(abs(self.obs.centroid.px.std_err_sq[1] - 2) < eps)
    assert(abs(self.obs.centroid.px.std_err_sq[2] - 3) < eps)

    assert(abs(self.obs.centroid.mm.position[0] - 100.1) < eps)
    assert(abs(self.obs.centroid.mm.position[1] - 200.1) < eps)
    assert(abs(self.obs.centroid.mm.position[2] - 300.1) < eps)
    assert(abs(self.obs.centroid.mm.variance[0] - 10.1) < eps)
    assert(abs(self.obs.centroid.mm.variance[1] - 20.1) < eps)
    assert(abs(self.obs.centroid.mm.variance[2] - 30.1) < eps)
    assert(abs(self.obs.centroid.mm.std_err_sq[0] - 1.1) < eps)
    assert(abs(self.obs.centroid.mm.std_err_sq[1] - 2.1) < eps)
    assert(abs(self.obs.centroid.mm.std_err_sq[2] - 3.1) < eps)

    assert(abs(self.obs.intensity.observed.value - 1234) < eps)
    assert(abs(self.obs.intensity.observed.variance - 42) < eps)
    assert(abs(self.obs.intensity.corrected.value - 1234.1) < eps)
    assert(abs(self.obs.intensity.corrected.variance - 42.1) < eps)


  def tst_equality(self):
    from dials.model.data import Observation

    obs2 = Observation(self.obs)
    assert(obs2 == self.obs)
    obs2.panel += 1
    assert(obs2 != self.obs)
    obs2 = Observation(self.obs)
    obs2.centroid.px.position = (0, 0, 0)
    assert(obs2 != self.obs)
    obs2 = Observation(self.obs)
    obs2.centroid.px.variance = (0, 0, 0)
    assert(obs2 != self.obs)
    obs2 = Observation(self.obs)
    obs2.centroid.px.std_err_sq = (0, 0, 0)
    assert(obs2 != self.obs)
    obs2 = Observation(self.obs)
    obs2.centroid.mm.position = (0, 0, 0)
    assert(obs2 != self.obs)
    obs2 = Observation(self.obs)
    obs2.centroid.mm.variance = (0, 0, 0)
    assert(obs2 != self.obs)
    obs2 = Observation(self.obs)
    obs2.centroid.mm.std_err_sq = (0, 0, 0)
    assert(obs2 != self.obs)
    obs2 = Observation(self.obs)
    obs2.intensity.observed.value = 0
    assert(obs2 != self.obs)
    obs2 = Observation(self.obs)
    obs2.intensity.observed.variance = 0
    assert(obs2 != self.obs)
    obs2 = Observation(self.obs)
    obs2.intensity.corrected.value = 0
    assert(obs2 != self.obs)
    obs2 = Observation(self.obs)
    obs2.intensity.corrected.variance = 0
    assert(obs2 != self.obs)

if __name__ == '__main__':
  from dials.test import cd_auto
  with cd_auto(__file__):
    test = Test()
    test.run()
