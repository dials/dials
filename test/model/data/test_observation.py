from __future__ import absolute_import, division, print_function

import pytest

@pytest.fixture
def observation():
  from dials.model.data import Observation
  obs = Observation()
  obs.panel = 1
  obs.centroid.px.position = (100, 200, 300)
  obs.centroid.px.variance = (10, 20, 30)
  obs.centroid.px.std_err_sq = (1, 2, 3)
  obs.centroid.mm.position = (100.1, 200.1, 300.1)
  obs.centroid.mm.variance = (10.1, 20.1, 30.1)
  obs.centroid.mm.std_err_sq = (1.1, 2.1, 3.1)
  obs.intensity.observed.value = 1234
  obs.intensity.observed.variance = 42
  obs.intensity.corrected.value = 1234.1
  obs.intensity.corrected.variance = 42.1
  return obs

def test_data(observation):
  eps = 1e-7

  assert observation.panel == 1

  assert abs(observation.centroid.px.position[0] - 100) < eps
  assert abs(observation.centroid.px.position[1] - 200) < eps
  assert abs(observation.centroid.px.position[2] - 300) < eps
  assert abs(observation.centroid.px.variance[0] - 10) < eps
  assert abs(observation.centroid.px.variance[1] - 20) < eps
  assert abs(observation.centroid.px.variance[2] - 30) < eps
  assert abs(observation.centroid.px.std_err_sq[0] - 1) < eps
  assert abs(observation.centroid.px.std_err_sq[1] - 2) < eps
  assert abs(observation.centroid.px.std_err_sq[2] - 3) < eps

  assert abs(observation.centroid.mm.position[0] - 100.1) < eps
  assert abs(observation.centroid.mm.position[1] - 200.1) < eps
  assert abs(observation.centroid.mm.position[2] - 300.1) < eps
  assert abs(observation.centroid.mm.variance[0] - 10.1) < eps
  assert abs(observation.centroid.mm.variance[1] - 20.1) < eps
  assert abs(observation.centroid.mm.variance[2] - 30.1) < eps
  assert abs(observation.centroid.mm.std_err_sq[0] - 1.1) < eps
  assert abs(observation.centroid.mm.std_err_sq[1] - 2.1) < eps
  assert abs(observation.centroid.mm.std_err_sq[2] - 3.1) < eps

  assert abs(observation.intensity.observed.value - 1234) < eps
  assert abs(observation.intensity.observed.variance - 42) < eps
  assert abs(observation.intensity.corrected.value - 1234.1) < eps
  assert abs(observation.intensity.corrected.variance - 42.1) < eps

def test_equality(observation):
  from dials.model.data import Observation

  obs2 = Observation(observation)
  assert obs2 == observation
  obs2.panel += 1
  assert obs2 != observation
  obs2 = Observation(observation)
  obs2.centroid.px.position = (0, 0, 0)
  assert obs2 != observation
  obs2 = Observation(observation)
  obs2.centroid.px.variance = (0, 0, 0)
  assert obs2 != observation
  obs2 = Observation(observation)
  obs2.centroid.px.std_err_sq = (0, 0, 0)
  assert obs2 != observation
  obs2 = Observation(observation)
  obs2.centroid.mm.position = (0, 0, 0)
  assert obs2 != observation
  obs2 = Observation(observation)
  obs2.centroid.mm.variance = (0, 0, 0)
  assert obs2 != observation
  obs2 = Observation(observation)
  obs2.centroid.mm.std_err_sq = (0, 0, 0)
  assert obs2 != observation
  obs2 = Observation(observation)
  obs2.intensity.observed.value = 0
  assert obs2 != observation
  obs2 = Observation(observation)
  obs2.intensity.observed.variance = 0
  assert obs2 != observation
  obs2 = Observation(observation)
  obs2.intensity.corrected.value = 0
  assert obs2 != observation
  obs2 = Observation(observation)
  obs2.intensity.corrected.variance = 0
  assert obs2 != observation
