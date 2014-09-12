
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
    self.tst_import_integrate_hkl()
    self.tst_import_spot_xds()
    self.tst_from_xds_files()

  def tst_import_integrate_hkl(self):
    from dials.array_family import flex # import dependency
    from os.path import join
    from libtbx import easy_run

    # Call dials.import_xds
    easy_run.fully_buffered([
      'dials.import_xds',
      'input.method=reflections',
      join(self.path, 'INTEGRATE.HKL'),
      join(self.path, "experiments.json")
    ]).raise_if_errors()

    import cPickle as pickle
    table = pickle.load(open('integrate_hkl.pickle', 'rb'))

    assert('miller_index' in table)
    assert('id' in table)
    assert('panel' in table)
    assert('xyzcal.px' in table)
    assert('xyzobs.px.value' in table)
    assert('intensity.cor.value' in table)
    assert('intensity.cor.variance' in table)
    assert(len(table) == 174911)
    print 'OK'

  def tst_import_spot_xds(self):
    from os.path import join
    from libtbx import easy_run

    # Call dials.import_xds
    easy_run.fully_buffered([
      'dials.import_xds',
      'input.method=reflections',
      join(self.path, 'SPOT.XDS'),
    ]).raise_if_errors()

    import cPickle as pickle
    table = pickle.load(open('spot_xds.pickle', 'rb'))

    assert('miller_index' in table)
    assert('id' in table)
    assert('panel' in table)
    assert('xyzobs.px.value' in table)
    assert('intensity.sum.value' in table)
    assert(len(table) == 742)
    print 'OK'

    # Call dials.import_xds
    easy_run.fully_buffered([
      'dials.import_xds',
      'input.method=reflections',
      join(self.path, 'SPOT.XDS'),
      'remove_invalid=True',
    ]).raise_if_errors()

    import cPickle as pickle
    table = pickle.load(open('spot_xds.pickle', 'rb'))

    assert('miller_index' in table)
    assert('id' in table)
    assert('panel' in table)
    assert('xyzobs.px.value' in table)
    assert('intensity.sum.value' in table)
    assert(len(table) == 664)
    print 'OK'


  def tst_from_xds_files(self):
    from subprocess import call
    from os.path import abspath, exists
    from os import chdir

    # Import from the image files
    path = abspath(self.path)
    chdir(path)
    call('dials.import_xds ./ output=import_experiments.json > /dev/null', shell=True)

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
  from dials.test import cd_auto
  with cd_auto(__file__):
    test = Test()
    test.run()
