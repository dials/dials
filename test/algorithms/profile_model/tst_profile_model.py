
from __future__ import absolute_import, division, print_function

class Test(object):

  def __init__(self):
    pass

  def run(self):
    self.tst_load_and_dump()

  def tst_load_and_dump(self):
    from dials.algorithms.profile_model.gaussian_rs import Model
    d1 = {
      '__id__' : 'gaussian_rs',
      'n_sigma' : 3,
      'sigma_b' : 1,
      'sigma_m' : 2
    }
    d2 = {
      '__id__' : 'gaussian_rs',
      'n_sigma' : 2,
      'sigma_b' : 4,
      'sigma_m' : 5
    }
    model1 = Model.from_dict(d1)
    model2 = Model.from_dict(d2)
    assert(model1.n_sigma() == 3)
    assert(model1.sigma_b() == 1)
    assert(model1.sigma_m() == 2)
    assert(model2.n_sigma() == 2)
    assert(model2.sigma_b() == 4)
    assert(model2.sigma_m() == 5)

if __name__ == '__main__':
  test = Test()
  test.run()
