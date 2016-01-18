from __future__ import division

# FIXME
# Currently looking at % difference between predicted background in foreground
# region. Should probably calculate the errors properly and look at z score
# distribution or something

class TestConstant2d(object):

  def __init__(self):
    from os.path import join, isfile
    import libtbx.load_env
    try:
      dials_regression = libtbx.env.dist_path('dials_regression')
    except KeyError, e:
      print 'SKIP: dials_regression not configured'
      exit(0)

    # The base path
    path = join(dials_regression, 'background_test_data', 'simulated')

    # The paths to the reflection files
    self.refl_filenames = [
      join(path, 'simulated_n1000_i0_b10.pickle'),
      join(path, 'simulated_n1000_i0_b100.pickle'),
      join(path, 'simulated_n1000_i0_b1000.pickle'),
      join(path, 'simulated_r_n1000_i0_b1000.pickle'),

      #join(path, 'simulated_n1000_i100_b0.pickle'),
      join(path, 'simulated_n1000_i100_b10.pickle'),
      join(path, 'simulated_n1000_i100_b100.pickle'),
      join(path, 'simulated_n1000_i100_b1000.pickle'),

      #join(path, 'simulated_n1000_i10000_b0.pickle'),
      join(path, 'simulated_n1000_i10000_b10.pickle'),
      join(path, 'simulated_n1000_i10000_b100.pickle'),
      join(path, 'simulated_n1000_i10000_b1000.pickle'),
    ]

    # Check the files exist
    for filename in self.refl_filenames:
      if not isfile(filename):
        print 'SKIP: simulated test data does not exist'
        print 'Generate by running the following commands:'
        print ' cd dials_regression/background_test_data/simulated'
        print ' ./simulate'
        exit(0)

  def run(self):
    from dials.algorithms.background.simple import Creator
    from dials.algorithms.background.simple import Constant2dModeller
    from dials.algorithms.background.simple import TruncatedOutlierRejector
    from dials.algorithms.background.simple import NSigmaOutlierRejector
    from dials.algorithms.background.simple import NormalOutlierRejector

    modeller = Constant2dModeller()

    outlier_rejector = [
      None,
      TruncatedOutlierRejector(0.01, 0.01),
      NSigmaOutlierRejector(3.0, 3.0),
      NormalOutlierRejector(10)
    ]

    for rejector in outlier_rejector:
      self.tst(Creator(modeller, rejector))

  def tst(self, creator):
    for filename in self.refl_filenames:
      self.tst_for_dataset(creator, filename)
    print 'OK'

  def tst_for_dataset(self, creator, filename):
    from dials.array_family import flex
    from dials.algorithms.shoebox import MaskCode
    print filename
    rlist = flex.reflection_table.from_pickle(filename)
    shoebox = rlist['shoebox']
    background = [sb.background.deep_copy() for sb in shoebox]
    success = creator(shoebox)
    assert(success.count(True) == len(success))
    diff = []
    for i in range(len(rlist)):
      mask = flex.bool([(m & MaskCode.Foreground) != 0 for m in shoebox[i].mask])
      px1 = background[i].select(mask)
      px2 = shoebox[i].background.select(mask)
      den = max([flex.mean(px1), 1.0])
      diff.append(flex.mean(px2 - px1) / den)
    diff = flex.double(diff)
    mv = flex.mean_and_variance(flex.double(diff))
    mean = mv.mean()
    sdev = mv.unweighted_sample_standard_deviation()
    try:
      assert(abs(mean) < 0.01)
    except Exception:
      print "Mean: %f, Sdev: %f", mean, sdev
      from matplotlib import pylab
      pylab.hist(diff)
      pylab.show()
      raise

class TestConstant3d(object):

  def __init__(self):
    from os.path import join, isfile
    import libtbx.load_env
    try:
      dials_regression = libtbx.env.dist_path('dials_regression')
    except KeyError, e:
      print 'SKIP: dials_regression not configured'
      exit(0)

    # The base path
    path = join(dials_regression, 'background_test_data', 'simulated')

    # The paths to the reflection files
    self.refl_filenames = [
      join(path, 'simulated_n1000_i0_b10.pickle'),
      join(path, 'simulated_n1000_i0_b100.pickle'),
      join(path, 'simulated_n1000_i0_b1000.pickle'),
      join(path, 'simulated_r_n1000_i0_b1000.pickle'),

      #join(path, 'simulated_n1000_i100_b0.pickle'),
      join(path, 'simulated_n1000_i100_b10.pickle'),
      join(path, 'simulated_n1000_i100_b100.pickle'),
      join(path, 'simulated_n1000_i100_b1000.pickle'),

      #join(path, 'simulated_n1000_i10000_b0.pickle'),
      join(path, 'simulated_n1000_i10000_b10.pickle'),
      join(path, 'simulated_n1000_i10000_b100.pickle'),
      join(path, 'simulated_n1000_i10000_b1000.pickle'),
    ]

    # Check the files exist
    for filename in self.refl_filenames:
      if not isfile(filename):
        print 'SKIP: simulated test data does not exist'
        print 'Generate by running the following commands:'
        print ' cd dials_regression/background_test_data/simulated'
        print ' ./simulate'
        exit(0)

  def run(self):
    from dials.algorithms.background.simple import Creator
    from dials.algorithms.background.simple import Constant3dModeller
    from dials.algorithms.background.simple import TruncatedOutlierRejector
    from dials.algorithms.background.simple import NSigmaOutlierRejector
    from dials.algorithms.background.simple import NormalOutlierRejector

    modeller = Constant3dModeller()

    outlier_rejector = [
      None,
      TruncatedOutlierRejector(0.01, 0.01),
      NSigmaOutlierRejector(3.0, 3.0),
      NormalOutlierRejector(10)
    ]

    for rejector in outlier_rejector:
      self.tst(Creator(modeller, rejector))

  def tst(self, creator):
    for filename in self.refl_filenames:
      self.tst_for_dataset(creator, filename)
    print 'OK'

  def tst_for_dataset(self, creator, filename):
    from dials.array_family import flex
    from dials.algorithms.shoebox import MaskCode
    print filename
    rlist = flex.reflection_table.from_pickle(filename)
    shoebox = rlist['shoebox']
    background = [sb.background.deep_copy() for sb in shoebox]
    success = creator(shoebox)
    assert(success.count(True) == len(success))
    diff = []
    for i in range(len(rlist)):
      mask = flex.bool([(m & MaskCode.Foreground) != 0 for m in shoebox[i].mask])
      px1 = background[i].select(mask)
      px2 = shoebox[i].background.select(mask)
      den = max([flex.mean(px1), 1.0])
      diff.append(flex.mean(px2 - px1) / den)
    diff = flex.double(diff)
    mv = flex.mean_and_variance(flex.double(diff))
    mean = mv.mean()
    sdev = mv.unweighted_sample_standard_deviation()
    try:
      assert(abs(mean) < 0.01)
    except Exception:
      print "Mean: %f, Sdev: %f", mean, sdev
      #from matplotlib import pylab
      #pylab.hist(diff)
      #pylab.show()
      raise


class TestLinear2d(object):

  def __init__(self):
    from os.path import join, isfile
    import libtbx.load_env
    try:
      dials_regression = libtbx.env.dist_path('dials_regression')
    except KeyError, e:
      print 'SKIP: dials_regression not configured'
      exit(0)

    # The base path
    path = join(dials_regression, 'background_test_data', 'simulated')

    # The paths to the reflection files
    self.refl_filenames = [
      join(path, 'simulated_n1000_i0_b10,1,1.pickle'),
      join(path, 'simulated_n1000_i0_b100,1,1.pickle'),
      join(path, 'simulated_n1000_i0_b1000,1,1.pickle'),
      join(path, 'simulated_r_n1000_i0_b1000,1,1.pickle'),

      #join(path, 'simulated_n1000_i100_b0,0,0.pickle'),
      join(path, 'simulated_n1000_i100_b10,1,1.pickle'),
      join(path, 'simulated_n1000_i100_b100,1,1.pickle'),
      join(path, 'simulated_n1000_i100_b1000,1,1.pickle'),

      #join(path, 'simulated_n1000_i10000_b0,0,0.pickle'),
      join(path, 'simulated_n1000_i10000_b10,1,1.pickle'),
      join(path, 'simulated_n1000_i10000_b100,1,1.pickle'),
      join(path, 'simulated_n1000_i10000_b1000,1,1.pickle'),
    ]

    # Check the files exist
    for filename in self.refl_filenames:
      if not isfile(filename):
        print 'SKIP: simulated test data does not exist'
        print 'Generate by running the following commands:'
        print ' cd dials_regression/background_test_data/simulated'
        print ' ./simulate'
        exit(0)

  def run(self):
    from dials.algorithms.background.simple import Creator
    from dials.algorithms.background.simple import Linear2dModeller
    from dials.algorithms.background.simple import TruncatedOutlierRejector
    from dials.algorithms.background.simple import NSigmaOutlierRejector
    from dials.algorithms.background.simple import NormalOutlierRejector

    modeller = Linear2dModeller()

    outlier_rejector = [
      None,
      TruncatedOutlierRejector(0.01, 0.01),
      NSigmaOutlierRejector(3.0, 3.0),
      NormalOutlierRejector(10)
    ]

    for rejector in outlier_rejector:
      self.tst(Creator(modeller, rejector))

  def tst(self, creator):
    for filename in self.refl_filenames:
      self.tst_for_dataset(creator, filename)
    print 'OK'

  def tst_for_dataset(self, creator, filename):
    from dials.array_family import flex
    from dials.algorithms.shoebox import MaskCode
    print filename
    rlist = flex.reflection_table.from_pickle(filename)
    shoebox = rlist['shoebox']
    background = [sb.background.deep_copy() for sb in shoebox]
    success = creator(shoebox)
    assert(success.count(True) == len(success))
    diff = []
    for i in range(len(rlist)):
      mask = flex.bool([(m & MaskCode.Foreground) != 0 for m in shoebox[i].mask])
      px1 = background[i].select(mask)
      px2 = shoebox[i].background.select(mask)
      den = max([flex.mean(px1), 1.0])
      diff.append(flex.mean(px2 - px1) / den)
    diff = flex.double(diff)
    mv = flex.mean_and_variance(flex.double(diff))
    mean = mv.mean()
    sdev = mv.unweighted_sample_standard_deviation()
    try:
      assert(abs(mean) < 0.01)
    except Exception:
      print "Mean: %f, Sdev: %f", mean, sdev
      #from matplotlib import pylab
      #pylab.hist(diff)
      #pylab.show()
      raise


class TestLinear3d(object):

  def __init__(self):
    from os.path import join, isfile
    import libtbx.load_env
    try:
      dials_regression = libtbx.env.dist_path('dials_regression')
    except KeyError, e:
      print 'SKIP: dials_regression not configured'
      exit(0)

    # The base path
    path = join(dials_regression, 'background_test_data', 'simulated')

    # The paths to the reflection files
    self.refl_filenames = [
      join(path, 'simulated_n1000_i0_b10,1,1,1.pickle'),
      join(path, 'simulated_n1000_i0_b100,1,1,1.pickle'),
      join(path, 'simulated_n1000_i0_b1000,1,1,1.pickle'),
      join(path, 'simulated_r_n1000_i0_b1000,1,1,1.pickle'),

      #join(path, 'simulated_n1000_i100_b0,0,0,0.pickle'),
      join(path, 'simulated_n1000_i100_b10,1,1,1.pickle'),
      join(path, 'simulated_n1000_i100_b100,1,1,1.pickle'),
      join(path, 'simulated_n1000_i100_b1000,1,1,1.pickle'),

      #join(path, 'simulated_n1000_i10000_b0,0,0,0.pickle'),
      join(path, 'simulated_n1000_i10000_b10,1,1,1.pickle'),
      join(path, 'simulated_n1000_i10000_b100,1,1,1.pickle'),
      join(path, 'simulated_n1000_i10000_b1000,1,1,1.pickle'),
    ]

    # Check the files exist
    for filename in self.refl_filenames:
      if not isfile(filename):
        print 'SKIP: simulated test data does not exist'
        print 'Generate by running the following commands:'
        print ' cd dials_regression/background_test_data/simulated'
        print ' ./simulate'
        exit(0)

  def run(self):
    from dials.algorithms.background.simple import Creator
    from dials.algorithms.background.simple import Linear3dModeller
    from dials.algorithms.background.simple import TruncatedOutlierRejector
    from dials.algorithms.background.simple import NSigmaOutlierRejector
    from dials.algorithms.background.simple import NormalOutlierRejector

    modeller = Linear3dModeller()

    outlier_rejector = [
      None,
      TruncatedOutlierRejector(0.01, 0.01),
      NSigmaOutlierRejector(3.0, 3.0),
      NormalOutlierRejector(10)
    ]

    for rejector in outlier_rejector:
      self.tst(Creator(modeller, rejector))

  def tst(self, creator):
    for filename in self.refl_filenames:
      self.tst_for_dataset(creator, filename)
    print 'OK'

  def tst_for_dataset(self, creator, filename):
    from dials.array_family import flex
    from dials.algorithms.shoebox import MaskCode
    print filename
    rlist = flex.reflection_table.from_pickle(filename)
    shoebox = rlist['shoebox']

    # FIXME doesn't work for single image
    zr = flex.int([s.bbox[5]-s.bbox[4] for s in shoebox])
    rlist = rlist.select(zr > 1)
    shoebox = rlist['shoebox']

    background = [sb.background.deep_copy() for sb in shoebox]
    success = creator(shoebox)
    assert(success.count(True) == len(success))
    diff = []
    for i in range(len(rlist)):
      mask = flex.bool([(m & MaskCode.Foreground) != 0 for m in shoebox[i].mask])
      px1 = background[i].select(mask)
      px2 = shoebox[i].background.select(mask)
      den = max([flex.mean(px1), 1.0])
      diff.append(flex.mean(px2 - px1) / den)
    diff = flex.double(diff)
    mv = flex.mean_and_variance(flex.double(diff))
    mean = mv.mean()
    sdev = mv.unweighted_sample_standard_deviation()
    try:
      assert(abs(mean) < 0.01)
    except Exception:
      print "Mean: %f, Sdev: %f", mean, sdev
      #from matplotlib import pylab
      #pylab.hist(diff)
      #pylab.show()
      raise


if __name__ == '__main__':

  test = TestConstant2d()
  test.run()

  test = TestConstant3d()
  test.run()

  test = TestLinear2d()
  test.run()

  test = TestLinear3d()
  test.run()
