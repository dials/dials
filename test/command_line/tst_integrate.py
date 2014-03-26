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

  def run(self):
    from os.path import abspath, join
    from libtbx import easy_run
    import os
    from uuid import uuid4

    dirname ='tmp_%s' % uuid4().hex
    os.mkdir(dirname)
    os.chdir(dirname)

    # Call dials.integrate
    easy_run.fully_buffered([
      'dials.integrate',
      join(self.path, 'experiments.json'),
      'integration.intensity.algorithm=sum3d',
      'integration.shoebox.sigma_b=0.058',
      'integration.shoebox.sigma_m=0.157',
    ]).raise_if_errors()

    import cPickle as pickle
    table = pickle.load(open('integrated.pickle', 'rb'))
    assert(len(table) == 757)
    assert('id' in table)
    for row in table:
      assert(row['id'] == 0)
    print 'OK'


if __name__ == '__main__':
  from dials.test import cd_auto
  with cd_auto(__file__):
    test = Test()
    test.run()
