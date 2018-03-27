from __future__ import absolute_import, division, print_function

class Test(object):

  def __init__(self):
    from dials.model.data import Prediction
    self.pred = Prediction()
    self.pred.miller_index = (1, 2, 3)
    self.pred.beam_vector = (4, 5, 6)
    self.pred.panel = 7
    self.pred.entering = True
    self.pred.position.px = (8, 9, 10)
    self.pred.position.mm = (11, 12, 13)

  def run(self):
    self.tst_data()
    self.tst_equality()

  def tst_data(self):
    eps = 1e-7

    assert(self.pred.miller_index == (1, 2, 3))
    assert(self.pred.panel == 7)
    assert(self.pred.entering)
    assert(abs(self.pred.beam_vector[0] - 4) <= eps)
    assert(abs(self.pred.beam_vector[1] - 5) <= eps)
    assert(abs(self.pred.beam_vector[2] - 6) <= eps)
    assert(abs(self.pred.position.px[0] - 8) <= eps)
    assert(abs(self.pred.position.px[1] - 9) <= eps)
    assert(abs(self.pred.position.px[2] - 10) <= eps)
    assert(abs(self.pred.position.mm[0] - 11) <= eps)
    assert(abs(self.pred.position.mm[1] - 12) <= eps)
    assert(abs(self.pred.position.mm[2] - 13) <= eps)


  def tst_equality(self):
    from dials.model.data import Prediction

    pred2 = Prediction(self.pred)
    assert(pred2 == self.pred)
    pred2 = Prediction(self.pred)
    pred2.miller_index = (0, 0, 0)
    assert(pred2 != self.pred)
    pred2 = Prediction(self.pred)
    pred2.beam_vector = (0, 0, 0)
    assert(pred2 != self.pred)
    pred2 = Prediction(self.pred)
    pred2.panel = 0
    assert(pred2 != self.pred)
    pred2 = Prediction(self.pred)
    pred2.entering = False
    assert(pred2 != self.pred)
    pred2 = Prediction(self.pred)
    pred2.position.px = (0, 0, 0)
    assert(pred2 != self.pred)
    pred2 = Prediction(self.pred)
    pred2.position.mm = (0, 0, 0)
    assert(pred2 != self.pred)

if __name__ == '__main__':
  from dials.test import cd_auto
  with cd_auto(__file__):
    test = Test()
    test.run()
