from __future__ import division


class Test(object):

  def __init__(self):
    from os.path import join
    from dials.array_family import flex
    import libtbx.load_env
    try:
      dials_regression = libtbx.env.dist_path('dials_regression')
    except KeyError, e:
      print 'FAIL: dials_regression not configured'
      exit(0)

    self.path = join(dials_regression, "centroid_test_data")
    #self.path = join(dials_regression, "integration_test_data")
  def run(self):
    from os.path import abspath, join
    from libtbx import easy_run

    # Call dials.integrate
    easy_run.fully_buffered([
      'dials.integrate',
      join(self.path, 'experiments.json'),
      'integration.algorithm=sum3d',
    ]).raise_if_errors()

    import cPickle as pickle
    table = pickle.load(open('integrated.pickle', 'rb'))
    assert(len(table) == 360)
    print 'OK'

if __name__ == '__main__':
  test = Test()
  test.run()
