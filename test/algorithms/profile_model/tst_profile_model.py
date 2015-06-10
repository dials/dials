
from __future__ import division

class Test(object):

  def __init__(self):
    pass

  def run(self):
    self.tst_load_and_dump()

  def tst_load_and_dump(self):
    from dials.algorithms.profile_model.gaussian_rs import Factory
    from dials.algorithms.profile_model.factory import phil_scope
    from libtbx.phil import parse

    user_phil = parse('''
      profile {
        gaussian_rs {
          model {
            sigma_b = 1
            sigma_m = 2
          }
          model {
            n_sigma = 2
            sigma_b = 4
            sigma_m = 5
          }
        }
      }
      ''')
    params = phil_scope.fetch(source=user_phil).extract()
    model = Factory.load(params.profile)
    assert(len(model) == 2)
    assert(model[0].n_sigma() == 3)
    assert(model[0].sigma_b() == 1)
    assert(model[0].sigma_m() == 2)
    assert(model[1].n_sigma() == 2)
    assert(model[1].sigma_b() == 4)
    assert(model[1].sigma_m() == 5)
    print 'OK'

    model_phil = model.dump().fetch_diff(user_phil)
    assert(model_phil.as_str() == '')
    print 'OK'

if __name__ == '__main__':
  test = Test()
  test.run()
