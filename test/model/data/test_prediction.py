from __future__ import absolute_import, division, print_function

import pytest

@pytest.fixture
def prediction():
  from dials.model.data import Prediction
  pred = Prediction()
  pred.miller_index = (1, 2, 3)
  pred.beam_vector = (4, 5, 6)
  pred.panel = 7
  pred.entering = True
  pred.position.px = (8, 9, 10)
  pred.position.mm = (11, 12, 13)
  return pred

def test_data(prediction):
  eps = 1e-7

  assert prediction.miller_index == (1, 2, 3)
  assert prediction.panel == 7
  assert prediction.entering
  assert abs(prediction.beam_vector[0] - 4) <= eps
  assert abs(prediction.beam_vector[1] - 5) <= eps
  assert abs(prediction.beam_vector[2] - 6) <= eps
  assert abs(prediction.position.px[0] - 8) <= eps
  assert abs(prediction.position.px[1] - 9) <= eps
  assert abs(prediction.position.px[2] - 10) <= eps
  assert abs(prediction.position.mm[0] - 11) <= eps
  assert abs(prediction.position.mm[1] - 12) <= eps
  assert abs(prediction.position.mm[2] - 13) <= eps

def test_equality(prediction):
  from dials.model.data import Prediction

  pred2 = Prediction(prediction)
  assert pred2 == prediction
  pred2 = Prediction(prediction)
  pred2.miller_index = (0, 0, 0)
  assert pred2 != prediction
  pred2 = Prediction(prediction)
  pred2.beam_vector = (0, 0, 0)
  assert pred2 != prediction
  pred2 = Prediction(prediction)
  pred2.panel = 0
  assert pred2 != prediction
  pred2 = Prediction(prediction)
  pred2.entering = False
  assert pred2 != prediction
  pred2 = Prediction(prediction)
  pred2.position.px = (0, 0, 0)
  assert pred2 != prediction
  pred2 = Prediction(prediction)
  pred2.position.mm = (0, 0, 0)
  assert pred2 != prediction
