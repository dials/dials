
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
    self.tst_from_xds_files()

  def tst_import_integrate_hkl(self):

    from os.path import abspath, join
    from libtbx import easy_run

    # Call dials.integrate
    easy_run.fully_buffered([
      'dials.import_xds',
      '-i', 'reflections',
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
      '-i', 'reflections',
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
      '-i', 'reflections',
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


  def tst_from_xds_files(self):
    from subprocess import call
    import difflib
    from os.path import join, abspath, exists
    from os import chdir

    # Import from the image files
    path = abspath(self.path)
    chdir(path)
    call('dials.import_xds ./ -o import_experiments.json > /dev/null', shell=True)

    assert(exists("import_experiments.json"))

    # Get the expected output
    #expected = self.expected_import_from_xds_files()

    ## Read the created file and do a diff
    #with open("experiments.json", "r") as infile:
      #lines_a = infile.read().splitlines()
      #lines_a = [l.strip() for l in lines_a if "\"template\"" not in l]
      #diff = list(difflib.context_diff(
        #lines_a,
        #[l.strip() for l in expected.splitlines()]))
      #n = len(diff)
      #for i, line in enumerate(diff):
        #print line
      #assert(n == 0)

    print 'OK'
if __name__ == '__main__':
  test = Test()
  test.run()
