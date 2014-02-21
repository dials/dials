

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
    from os.path import abspath, join, exists
    from libtbx import easy_run

    # Call dials.extract
    easy_run.fully_buffered([
      'dials.extract',
      join(self.path, 'experiments.json'),
    ]).raise_if_errors()

    assert(exists("extracted.tar"))

    print 'OK'

if __name__ == '__main__':
  test = Test()
  test.run()
