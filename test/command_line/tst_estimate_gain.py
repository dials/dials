
from __future__ import division
from dials.array_family import flex # import dependency


class Test(object):

  def __init__(self):
    from os.path import join
    import libtbx.load_env
    try:
      dials_regression = libtbx.env.dist_path('dials_regression')
    except KeyError, e:
      print 'SKIP: dials_regression not configured'
      exit(0)

    self.path = join(dials_regression, "centroid_test_data")

  def run(self):

    from os.path import join, exists
    from libtbx import easy_run
    import os

    input_filename = join(self.path, "datablock.json")

    easy_run.fully_buffered(
      ['dials.estimate_gain',
       'input.datablock=%s' % input_filename]).raise_if_errors()

    print 'OK'

if __name__ == '__main__':
  from dials.test import cd_auto
  with cd_auto(__file__):
    test = Test()
    test.run()
