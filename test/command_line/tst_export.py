
from __future__ import division


class Test(object):

  def __init__(self):
    from os.path import join
    from libtbx import easy_run
    import libtbx.load_env
    try:
      dials_regression = libtbx.env.dist_path('dials_regression')
    except KeyError, e:
      print 'FAIL: dials_regression not configured'
      exit(0)

    try:
      import h5py # implicit import
    except ImportError:
      print "Skipping: can't import module h5py"
      exit(0)

    self.path = join(dials_regression, "centroid_test_data")

    # Call dials.integrate
    easy_run.fully_buffered([
      'dials.integrate',
      join(self.path, 'experiments.json'),
      join(self.path, 'profile.phil'),
      'profile.fitting=False',
    ]).raise_if_errors()

  def run(self):
    from libtbx import easy_run
    from dials.array_family import flex
    from dials.model.data import Shoebox

    # Call dials.export
    easy_run.fully_buffered([
      'dials.export',
      'integrated.pickle',
      'integrated.h5'
    ]).raise_if_errors()

    table1 = flex.reflection_table.from_pickle('integrated.pickle')
    table2 = flex.reflection_table.from_h5('integrated.h5')
    eps = 1e-7
    assert(table1.nrows() == table2.nrows())
    assert(table1.ncols() == table2.ncols())
    assert(sorted(table1.keys()) == sorted(table2.keys()))
    for key in table1.keys():
      col1 = table1[key]
      col2 = table2[key]
      assert(type(col1) == type(col1))
      for e1, e2 in zip(col1, col2):
        if isinstance(e1, tuple):
          assert(len(e1) == len(e2))
          for ee1, ee2 in zip(e1, e2):
            assert(abs(ee1 - ee2) < eps)
        elif isinstance(e1, Shoebox):
          pass
        else:
          assert(abs(e1 - e2) < eps)

    print 'OK'

if __name__ == '__main__':
  from dials.test import cd_auto
  with cd_auto(__file__):
    test = Test()
    test.run()
