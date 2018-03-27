
from __future__ import absolute_import, division

class Test(object):

  def __init__(self):
    from os.path import join
    import libtbx.load_env
    try:
      dials_regression = libtbx.env.dist_path('dials_regression')
    except KeyError:
      print 'FAIL: dials_regression not configured'
      exit(0)

    self.path = join(dials_regression, "centroid_test_data")

  def run(self):
    from os.path import join, exists
    from libtbx import easy_run
    from dials.array_family import flex


    assert(exists(join(self.path, "integrated.pickle")))

    input_filename = join(self.path, "integrated.pickle")

    easy_run.fully_buffered([
      'dev.dials.sort_reflections',
      input_filename,
      'key=intensity.sum.value',
      'output=sorted1.pickle',
    ]).raise_if_errors()

    assert(exists("sorted1.pickle"))

    sorted1 = flex.reflection_table.from_pickle("sorted1.pickle")
    self.assert_sorted(sorted1['intensity.sum.value'])

    easy_run.fully_buffered([
      'dev.dials.sort_reflections',
      input_filename,
      'output=sorted2.pickle',
      'key=intensity.sum.value',
      'reverse=True'
    ]).raise_if_errors()

    assert(exists("sorted2.pickle"))

    sorted1 = flex.reflection_table.from_pickle("sorted2.pickle")
    self.assert_sorted(sorted1['intensity.sum.value'], reverse=True)

    # test default sort on miller_index
    easy_run.fully_buffered([
      'dev.dials.sort_reflections',
      input_filename,
      'output=sorted3.pickle'
    ]).raise_if_errors()

    assert(exists("sorted3.pickle"))

    sorted1 = flex.reflection_table.from_pickle("sorted3.pickle")
    mi1 = sorted1['miller_index']
    orig = flex.reflection_table.from_pickle(input_filename)
    mi2 = flex.miller_index(sorted(orig['miller_index']))
    assert mi1.all_eq(mi2)


  def assert_sorted(self, x, reverse=False):
    assert(len(x) > 0)
    x0 = x[0]
    for i in range(1, len(x)):
      x1 = x[i]
      if reverse is True:
        assert(x1 <= x0)
      else:
        assert(x1 >= x0)
      x0 = x1


if __name__ == '__main__':
  from dials.test import cd_auto
  with cd_auto(__file__):
    test = Test()
    test.run()
