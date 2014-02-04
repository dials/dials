
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
    self.tst_import_integrate_hkl()
    self.tst_import_spot_xds()

  def tst_import_integrate_hkl(self):

    from os.path import abspath, join
    from libtbx import easy_run

    # Call dials.integrate
    easy_run.fully_buffered([
      'dials.import_xds',
      join(self.path, 'INTEGRATE.HKL'),
    ]).raise_if_errors()

    import cPickle as pickle
    table = pickle.load(open('integrate_hkl.pickle', 'rb'))

    assert('hkl' in table)
    assert('id' in table)
    assert('panel' in table)
    assert('xyzcal.px' in table)
    assert('xyzobs.px.value' in table)
    assert('intensity.cor.value' in table)
    assert('intensity.cor.variance' in table)
    assert(len(table) == 174911)
    print 'OK'

  def tst_import_spot_xds(self):
    from os.path import abspath, join
    from libtbx import easy_run

    # Call dials.integrate
    easy_run.fully_buffered([
      'dials.import_xds',
      join(self.path, 'SPOT.XDS'),
    ]).raise_if_errors()

    import cPickle as pickle
    table = pickle.load(open('spot_xds.pickle', 'rb'))

    assert('hkl' in table)
    assert('id' in table)
    assert('panel' in table)
    assert('xyzobs.px.value' in table)
    assert('intensity.raw.value' in table)
    assert(len(table) == 742)
    print 'OK'

    # Call dials.integrate
    easy_run.fully_buffered([
      'dials.import_xds',
      join(self.path, 'SPOT.XDS'),
      '-r',
    ]).raise_if_errors()

    import cPickle as pickle
    table = pickle.load(open('spot_xds.pickle', 'rb'))

    assert('hkl' in table)
    assert('id' in table)
    assert('panel' in table)
    assert('xyzobs.px.value' in table)
    assert('intensity.raw.value' in table)
    assert(len(table) == 664)
    print 'OK'

if __name__ == '__main__':
  test = Test()
  test.run()
