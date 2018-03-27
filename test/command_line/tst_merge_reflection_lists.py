
from __future__ import absolute_import, division

class Test(object):

  def __init__(self):
    from os.path import join
    from dials.array_family import flex
    import libtbx.load_env
    from dials.array_family import flex
    try:
      dials_regression = libtbx.env.dist_path('dials_regression')
    except KeyError:
      print 'FAIL: dials_regression not configured'
      exit(0)

    self.path = join(dials_regression, "centroid_test_data")

    table = flex.reflection_table()
    table['hkl'] = flex.miller_index(360)
    table['id'] = flex.int(360)
    table['intensity.sum.value'] = flex.double(360)
    table.as_pickle("temp1.pickle")
    table.as_pickle("temp2.pickle")

  def run(self):
    from libtbx import easy_run
    from dials.array_family import flex

    # Call dials.merge_reflection_lists
    easy_run.fully_buffered([
      'dev.dials.merge_reflection_lists',
      'temp1.pickle',
      'temp2.pickle',
      'method=update'
    ]).raise_if_errors()

    table = flex.reflection_table.from_pickle('merged.pickle')
    assert(len(table) == 360)

    # Call dials.merge_reflection_lists
    easy_run.fully_buffered([
      'dev.dials.merge_reflection_lists',
      'temp1.pickle',
      'temp2.pickle',
      'method=extend'
    ]).raise_if_errors()

    table = flex.reflection_table.from_pickle('merged.pickle')
    assert(len(table) == 720)


if __name__ == '__main__':
  from dials.test import cd_auto
  with cd_auto(__file__):
    test = Test()
    test.run()
