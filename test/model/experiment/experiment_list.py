

from __future__ import division

from dials.model.experiment.experiment_list import \
  ExperimentList, Experiment, ExperimentListFactory

class TestExperiment(object):

  def __init__(self, path):
    self.path = path

  def run(self):
    self.tst_contains()
    self.tst_equality()
    self.tst_consistent()

  def tst_contains(self):
    from dxtbx.model import Beam, Detector, Goniometer, Scan
    from dials.model.experiment import Crystal

    # Create a load of models
    b1 = Beam()
    d1 = Detector()
    g1 = Goniometer()
    s1 = Scan()
    c1 = Crystal((1, 0, 0), (0, 1, 0), (0, 0, 1), 0)

    # Create an experiment
    e = Experiment(
      beam=b1, detector=d1, goniometer=g1,
      scan=s1, crystal=c1, imageset=None)

    # Check experiment contains model
    assert(b1 in e)
    assert(d1 in e)
    assert(g1 in e)
    assert(s1 in e)
    assert(c1 in e)

    # Create a load of models that look the same but aren't
    b2 = Beam()
    d2 = Detector()
    g2 = Goniometer()
    s2 = Scan()
    c2 = Crystal((1, 0, 0), (0, 1, 0), (0, 0, 1), 0)

    # Check experiment doesn't contain model
    assert(b2 not in e)
    assert(d2 not in e)
    assert(g2 not in e)
    assert(s2 not in e)
    assert(c2 not in e)

    # Test passed
    print 'OK'

  def tst_equality(self):

    from dxtbx.model import Beam, Detector, Goniometer, Scan
    from dials.model.experiment import Crystal

    # Create a load of models
    b1 = Beam()
    d1 = Detector()
    g1 = Goniometer()
    s1 = Scan()
    c1 = Crystal((1, 0, 0), (0, 1, 0), (0, 0, 1), 0)

    # Create a load of models that look the same but aren't
    b2 = Beam()
    d2 = Detector()
    g2 = Goniometer()
    s2 = Scan()
    c2 = Crystal((1, 0, 0), (0, 1, 0), (0, 0, 1), 0)

    # Create an experiment
    e1 = Experiment(
      beam=b1, detector=d1, goniometer=g1,
      scan=s1, crystal=c1, imageset=None)

    # Create an experiment
    e2 = Experiment(
      beam=b1, detector=d1, goniometer=g1,
      scan=s1, crystal=c1, imageset=None)

    # Create an experiment
    e3 = Experiment(
      beam=b2, detector=d2, goniometer=g2,
      scan=s2, crystal=c2, imageset=None)

    # Check e1 equals e2 but not e3
    assert(e1 == e2)
    assert(e1 != e3)
    assert(e2 != e3)

    # Test passed
    print 'OK'

  def tst_consistent(self):

    from dxtbx.imageset2 import ImageSetFactory
    from glob import glob
    from os.path import join
    from dxtbx.model import Scan

    # Create a sweep
    sweep_filenames = join(self.path, 'centroid_test_data', 'centroid*.cbf')
    sweep = ImageSetFactory.new(glob(sweep_filenames))[0]

    # Create experiment with sweep and good scan
    e = Experiment(imageset=sweep, scan=sweep.get_scan())
    assert(e.is_consistent())

    # Create experiment with sweep and defective scan
    scan = sweep.get_scan()
    scan.set_image_range((1, 1))
    e = Experiment(imageset=sweep, scan=scan)
    assert(not e.is_consistent())

    ## Create experiment with imageset and good scan
    #assert(e.is_consistent())

    ## Create experiment with imageset and non-still scan
    #assert(not e.is_consistent())

    ## Create experiment with imageset and scan with more than 1 image
    #assert(not e.is_consistent())

    ## Create experiment with imageset and defective scan
    #assert(not e.is_consistent())

    # Test passed
    print 'OK'

class TestExperimentList(object):

  def __init__(self, path):
    pass

  def run(self):
    self.tst_contains()
    self.tst_index()
    self.tst_replace()
    self.tst_indices()
    self.tst_models()
    self.tst_to_dict()

  def tst_contains(self):
    print 'OK'

  def tst_index(self):
    print 'OK'

  def tst_replace(self):
    print 'OK'

  def tst_indices(self):
    print 'OK'

  def tst_models(self):
    print 'OK'

  def tst_to_dict(self):
    print 'OK'

class TestExperimentListFactory(object):

  def __init__(self, path):
    pass

  def run(self):
    self.tst_from_datablock()
    self.tst_from_json()
    self.tst_from_pickle()

  def tst_from_datablock(self):
    print 'OK'

  def tst_from_json(self):
    print 'OK'

  def tst_from_pickle(self):
    print 'OK'

class Test(object):
  def __init__(self):
    import libtbx

    dials_regression = libtbx.env.dist_path('dials_regression')
    if not dials_regression:
      print 'Skipping: dials_regresson not configured'
      exit(0)

    self.tst_experiment = TestExperiment(dials_regression)
    self.tst_list = TestExperimentList(dials_regression)
    self.tst_factory = TestExperimentListFactory(dials_regression)

  def run(self):
    self.tst_experiment.run()
    self.tst_list.run()
    self.tst_factory.run()

if __name__ == '__main__':
  test = Test()
  test.run()
