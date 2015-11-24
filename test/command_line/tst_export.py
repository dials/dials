
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
    self.experiments = join(self.path, "experiments.json")
    self.reflections = join(self.path, "integrated.pickle")

  def run(self):
    self.test_mtz()
    self.test_nxs()

  def test_nxs(self):
    from libtbx import easy_run
    from os.path import exists

    # Call dials.export
    easy_run.fully_buffered([
      'dials.export',
      'format=nxs',
      self.experiments,
      self.reflections
    ]).raise_if_errors()

    assert exists("hklout.nxs")

    print 'OK'

  def test_mtz(self):
    from libtbx import easy_run
    from os.path import exists

    # Call dials.export
    easy_run.fully_buffered([
      'dials.export',
      'format=mtz',
      self.experiments,
      self.reflections
    ]).raise_if_errors()

    assert exists("hklout.mtz")

    print 'OK'

if __name__ == '__main__':
  from dials.test import cd_auto
  with cd_auto(__file__):
    test = Test()
    test.run()
