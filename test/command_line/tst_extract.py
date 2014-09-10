

from __future__ import division


class Test(object):

  def __init__(self):
    from os.path import join
    import libtbx.load_env
    try:
      dials_regression = libtbx.env.dist_path('dials_regression')
    except KeyError, e:
      print 'FAIL: dials_regression not configured'
      exit(0)

    self.path = join(dials_regression, "centroid_test_data")

  def run(self):
    from os.path import join, exists
    from libtbx import easy_run

    # Call dials.extract
    easy_run.fully_buffered([
      'dials.extract',
      join(self.path, 'experiments.json'),
      join(self.path, 'profile.phil'),
    ]).raise_if_errors()

    assert(exists("shoebox.dat"))

    print 'OK'

if __name__ == '__main__':
  from dials.test import cd_auto
  with cd_auto(__file__):
    test = Test()
    test.run()
