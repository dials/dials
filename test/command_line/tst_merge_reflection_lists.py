
from __future__ import division


class Test(object):

  def __init__(self):
    from os.path import join
    from dials.array_family import flex
    import libtbx.load_env
    from libtbx import easy_run
    try:
      dials_regression = libtbx.env.dist_path('dials_regression')
    except KeyError, e:
      print 'FAIL: dials_regression not configured'
      exit(0)

    self.path = join(dials_regression, "centroid_test_data")

    # Call dials.integrate
    easy_run.fully_buffered([
      'dials.integrate',
      join(self.path, 'experiments.json'),
      'integration.algorithm=sum3d',
    ]).raise_if_errors()

  def run(self):
    from os.path import abspath, join
    from libtbx import easy_run
    from dials.array_family import flex

    # Call dials.integrate
    easy_run.fully_buffered([
      'dials.merge_reflection_lists',
      'integrated.pickle',
      'integrated.pickle',
      '-m', 'update'
    ]).raise_if_errors()

    table = flex.reflection_table.from_pickle('merged.pickle')
    assert(len(table) == 360)
    print 'OK'

    # Call dials.integrate
    easy_run.fully_buffered([
      'dials.merge_reflection_lists',
      'integrated.pickle',
      'integrated.pickle',
      '-m', 'extend'
    ]).raise_if_errors()

    table = flex.reflection_table.from_pickle('merged.pickle')
    assert(len(table) == 720)
    print 'OK'


if __name__ == '__main__':
  test = Test()
  test.run()
