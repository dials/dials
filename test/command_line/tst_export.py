
from __future__ import absolute_import, division


class Test(object):

  def __init__(self):
    from os.path import join
    from libtbx import easy_run
    import libtbx.load_env
    try:
      dials_regression = libtbx.env.dist_path('dials_regression')
    except KeyError:
      print 'FAIL: dials_regression not configured'
      exit(0)

    try:
      import h5py # implicit import
    except ImportError:
      print "Skipping: can't import module h5py"
      exit(0)

    self.path = join(dials_regression, "centroid_test_data")
    self.datablock = join(self.path, "datablock.json")
    self.strong = join(self.path, "strong.pickle")
    self.experiments = join(self.path, "experiments.json")
    self.reflections = join(self.path, "integrated.pickle")

  def run(self):
    self.test_mtz()
    self.test_nxs()
    self.test_xds_ascii()
    self.test_sadabs()
    self.test_json()

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

    assert exists("integrated.nxs")

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

    assert exists("integrated.mtz")

    print 'OK'

  def test_xds_ascii(self):
    from libtbx import easy_run
    from os.path import exists
    from libtbx.test_utils import approx_equal

    # Call dials.export
    easy_run.fully_buffered([
      'dials.export',
      'summation=true',
      'format=xds_ascii',
      self.experiments,
      self.reflections
    ]).raise_if_errors()

    assert exists("DIALS.HKL")

    psi_values = {
      (-9, 7, -10):153.430361,
      (-5, 11, -26):175.559441,
      (-4, 23, 24):129.468070,
      (2, 10, 20):147.947274
      }

    for record in open('DIALS.HKL', 'r'):
      if record.startswith('!'):
        continue
      tokens = record.split()
      hkl = tuple(map(int, tokens[:3]))
      if not hkl in psi_values:
        continue
      psi = float(tokens[-1])
      assert approx_equal(psi, psi_values[hkl], eps=0.1)

    print 'OK'

  def test_sadabs(self):
    from libtbx import easy_run
    from os.path import exists
    from libtbx.test_utils import approx_equal

    # Call dials.export
    easy_run.fully_buffered([
      'dials.export',
      'summation=true',
      'format=sadabs',
      self.experiments,
      self.reflections
    ]).raise_if_errors()

    assert exists("integrated.sad")

    direction_cosines = {
      (-9, 7, -10):(0.51253, -0.72107, 0.84696, -0.68476, -0.14130, -0.10561),
      (-5, 11, -26):(0.51310, -0.62895, 0.84711, -0.59223, -0.13830, -0.50366),
      (-4, 23, 24):(0.51308, -0.60578, 0.84711, -0.31416, -0.13840, 0.73099),
      (2, 10, 20):(0.51239, -0.46605, 0.84693, -0.61521, -0.14204, 0.63586)
      }

    for record in open('integrated.sad', 'r'):
      record = record.replace('-', ' -')
      tokens = record.split()
      hkl = tuple(map(int, tokens[:3]))
      cosines = tuple(map(float, tokens[6:12]))
      if not hkl in direction_cosines:
        continue
      assert approx_equal(cosines, direction_cosines[hkl], eps=0.001)

    print 'OK'

  def test_json(self):
    import json
    from libtbx import easy_run
    from os.path import exists
    from libtbx.test_utils import approx_equal
    from dxtbx.datablock import DataBlockFactory

    # Call dials.export
    easy_run.fully_buffered([
      'dials.export',
      'format=json',
      self.datablock,
      self.strong
    ]).raise_if_errors()


    assert exists('rlp.json')
    with open('rlp.json', 'rb') as f:
      d = json.load(f)
      assert d.keys() == ['imageset_id', 'datablocks', 'rlp', 'experiment_id'], d.keys()
      assert d['rlp'][:3] == [0.123454, 0.57687, 0.186465], d['rlp'][:3]
      assert d['imageset_id'][0] == 0
      assert d['experiment_id'][0] == 0
      assert len(d['datablocks']) == 1
      db = DataBlockFactory.from_dict(d['datablocks'])
      imgset = db[0].extract_imagesets()
      assert len(imgset) == 1

    # Call dials.export
    easy_run.fully_buffered([
      'dials.export',
      'format=json',
      self.experiments,
      self.reflections,
      'json.filename=integrated.json',
      'n_digits=4',
      'compact=False',
    ]).raise_if_errors()

    assert exists('integrated.json')
    with open('integrated.json', 'rb') as f:
      d = json.load(f)
      assert d.keys() == ['imageset_id', 'rlp', 'experiment_id'], d.keys()
      assert d['rlp'][:3] == [-0.5975, -0.6141, 0.4702], d['rlp'][:3]
      assert d['imageset_id'][0] == 0
      assert d['experiment_id'][0] == 0

    print 'OK'

if __name__ == '__main__':
  from dials.test import cd_auto
  with cd_auto(__file__):
    test = Test()
    test.run()
