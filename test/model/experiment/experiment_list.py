

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
    self.path = path
    self.el = self.generate()

  def run(self):
    self.tst_contains()
    self.tst_index()
    self.tst_replace()
    self.tst_indices()
    self.tst_models()
    self.tst_to_dict()

  def tst_contains(self):

    from dxtbx.model import Beam, Detector, Goniometer, Scan

    # Check all the models are found
    for e in self.el:
      assert(e.beam in self.el)
      assert(e.detector in self.el)
      assert(e.goniometer in self.el)
      assert(e.scan in self.el)

    # Create some more models
    b = Beam()
    d = Detector()
    g = Goniometer()
    s = Scan()

    # Check that models not in are not found
    assert(b not in self.el)
    assert(d not in self.el)
    assert(g not in self.el)
    assert(s not in self.el)

    # Test passed
    print 'OK'

  def tst_index(self):

    # Check the indices of exisiting experiments
    assert(self.el.index(self.el[0]) is 0)
    assert(self.el.index(self.el[1]) is 1)
    assert(self.el.index(self.el[2]) is 2)
    assert(self.el.index(self.el[3]) is 1)
    assert(self.el.index(self.el[4]) is 0)

    # Check index of non exisiting experiment
    try:
      self.el.index(Experiment())
      assert(False)
    except ValueError:
      pass

    # Test passed
    print 'OK'

  def tst_replace(self):

    # Get the models
    b = [e.beam for e in self.el]
    d = [e.detector for e in self.el]
    g = [e.goniometer for e in self.el]
    s = [e.scan for e in self.el]

    # Replace some models
    self.el.replace(b[0], b[1])
    assert(self.el[0].beam is b[1])
    assert(self.el[4].beam is b[1])

    # Replace again
    self.el[0].beam = b[0]
    self.el[4].beam = b[4]

    # Test passed
    print 'OK'

  def tst_indices(self):
    from dxtbx.model import Beam, Detector, Goniometer, Scan

    # Get the models
    b = [e.beam for e in self.el]
    d = [e.detector for e in self.el]
    g = [e.goniometer for e in self.el]
    s = [e.scan for e in self.el]

    # Check indices of beams
    assert(self.el.indices(b[0]) == [0, 4])
    assert(self.el.indices(b[1]) == [1, 3])
    assert(self.el.indices(b[2]) == [2])
    assert(self.el.indices(b[3]) == [1, 3])
    assert(self.el.indices(b[4]) == [0, 4])

    # Check indices of detectors
    assert(self.el.indices(d[0]) == [0, 4])
    assert(self.el.indices(d[1]) == [1, 3])
    assert(self.el.indices(d[2]) == [2])
    assert(self.el.indices(d[3]) == [1, 3])
    assert(self.el.indices(d[4]) == [0, 4])

    # Check indices of goniometer
    assert(self.el.indices(g[0]) == [0, 4])
    assert(self.el.indices(g[1]) == [1, 3])
    assert(self.el.indices(g[2]) == [2])
    assert(self.el.indices(g[3]) == [1, 3])
    assert(self.el.indices(g[4]) == [0, 4])

    # Check indices of scans
    assert(self.el.indices(s[0]) == [0, 4])
    assert(self.el.indices(s[1]) == [1, 3])
    assert(self.el.indices(s[2]) == [2])
    assert(self.el.indices(s[3]) == [1, 3])
    assert(self.el.indices(s[4]) == [0, 4])

    # Check some models not in the list
    assert(len(self.el.indices(Beam())) == 0)
    assert(len(self.el.indices(Detector())) == 0)
    assert(len(self.el.indices(Goniometer())) == 0)
    assert(len(self.el.indices(Scan())) == 0)

    # Test passed
    print 'OK'

  def tst_models(self):

    # Get all the unique models
    b = self.el.beams()
    d = self.el.detectors()
    g = self.el.goniometers()
    s = self.el.scans()

    # Check we have the expected number
    assert(len(b) == 3)
    assert(len(d) == 3)
    assert(len(g) == 3)
    assert(len(s) == 3)

    # Check we have the expected order
    assert(b[0] == self.el[0].beam)
    assert(b[1] == self.el[1].beam)
    assert(b[2] == self.el[2].beam)

    assert(d[0] == self.el[0].detector)
    assert(d[1] == self.el[1].detector)
    assert(d[2] == self.el[2].detector)

    assert(g[0] == self.el[0].goniometer)
    assert(g[0] == self.el[0].goniometer)
    assert(g[1] == self.el[1].goniometer)

    assert(s[2] == self.el[2].scan)
    assert(s[1] == self.el[1].scan)
    assert(s[2] == self.el[2].scan)

    # Test passed
    print 'OK'

  def tst_to_dict(self):

    # Convert the list to a dictionary
    obj = self.el.to_dict()

    # Check this is the right object
    assert(obj['__id__'] == 'ExperimentList')

    # Check length of items
    assert(len(obj['experiment']) == 5)
    assert(len(obj['beam']) == 3)
    assert(len(obj['detector']) == 3)
    assert(len(obj['goniometer']) == 3)
    assert(len(obj['scan']) == 3)

    # The expected models
    b = [0, 1, 2, 1, 0]
    d = [0, 1, 2, 1, 0]
    g = [0, 1, 2, 1, 0]
    s = [0, 1, 2, 1, 0]

    # Check all the experiments
    for i, eobj in enumerate(obj['experiment']):
      assert(eobj['__id__'] == 'Experiment')
      assert(eobj['beam'] == b[i])
      assert(eobj['detector'] == d[i])
      assert(eobj['goniometer'] == g[i])
      assert(eobj['scan'] == s[i])

    # Test passed
    print 'OK'

  def generate(self):
    from dxtbx.model import Beam, Detector, Goniometer, Scan

    # Initialise a list of experiments
    experiments = ExperimentList()

    # Create a few beams
    b1 = Beam()
    b2 = Beam()
    b3 = Beam()

    # Create a few detectors
    d1 = Detector()
    d2 = Detector()
    d3 = Detector()

    # Create a few goniometers
    g1 = Goniometer()
    g2 = Goniometer()
    g3 = Goniometer()

    # Create a few scans
    s1 = Scan()
    s2 = Scan()
    s3 = Scan()

    # Create a list of models
    b = [b1, b2, b3, b2, b1]
    d = [d1, d2, d3, d2, d1]
    g = [g1, g2, g3, g2, g1]
    s = [s1, s2, s3, s2, s1]

    # Populate with various experiments
    for i in range(5):
      experiments.append(Experiment(
        beam=b[i],
        detector=d[i],
        goniometer=g[i],
        scan=s[i]))

    # Return the list of experiments
    return experiments

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
