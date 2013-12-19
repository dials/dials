class Experiment(object):

  __slots__ = ('imageset', 'beam', 'detector', 'goniometer', 'scan', 'crystal')

  def __init__(self, imageset=None, beam=None, detector=None,
               goniometer=None, scan=None, crystal=None):
    self.imageset = imageset
    self.beam = beam
    self.detector = detector
    self.goniometer = goniometer
    self.scan = scan
    self.crystal = crystal

  def __contains__(self, item):
    return (item is self.imageset or
            item is self.beam or
            item is self.detector or
            item is self.goniometer or
            item is self.scan or
            item is self.crystal)

  def __eq__(self, other):
    if not isinstance(other, Experiment):
      return False
    return (self.imageset is other.imageset and
            self.beam is other.beam and
            self.detector is other.detector and
            self.goniometer is other.goniometer and
            self.scan is other.scan and
            self.crystal is other.crystal)

  def __ne__(self, other):
    return not self.__eq__(other)


class ExperimentManager(object):

  def __init__(self, item=None):
    if item is not None:
      self._data = list(item)
    else:
      self._data = list()

  def __setitem__(self, index, item):
    if isinstance(item, Experiment):
      self._data[index] = item
    else:
      raise TypeError('expected type Experiment, got %s' % type(item))

  def __getitem__(self, index):
    if isinstance(index, slice):
      return ExperimentManager(self._data[index])
    return self._data[index]

  def __delitem__(self, index):
    del self._data[index]

  def __len__(self):
    return len(self._data)

  def __iter__(self):
    for e in self._data:
      yield e

  def __contains__(self, item):
    return item in self._data or any(item in e for e in self._data)

  def index(self, item):
    return self._data.index(item)

  def append(self, item):
    if isinstance(item, Experiment):
      self._data.append(item)
    else:
      raise TypeError('expected type Experiment, got %s' % type(item))

  def extend(self, other):
    if isinstance(other, ExperimentManager):
      self._data.extend(other._data)
    else:
      raise TypeError('expected type ExperimentManager, got %s' % type(item))

  def replace(self, a, b):
    for i in self.indices(a):
      exp = self._data[i]
      if   exp.imageset is a:   exp.imageset = b
      elif exp.beam is a:       exp.beam = b
      elif exp.detector is a:   exp.detector = b
      elif exp.goniometer is a: exp.goniometer = b
      elif exp.scan is a:       exp.scan = b
      elif exp.crystal is a:    exp.crystal = b
      else: raise ValueError('unidentified model %s' % a)

  def remove(self, model):
    self.replace(model, None)

  def indices(self, model):
    if isinstance(model, list) or isinstance(model, tuple):
      return list(set.intersection(*[set(self.indices(m)) for m in model]))
    else:
      return [i for i, e in enumerate(self) if model in e]

  def beams(self):
    return list(set([e.beam for e in self]))

  def detectors(self):
    return list(set([e.detector for e in self]))

  def goniometers(self):
    return list(set([e.goniometer for e in self]))

  def scans(self):
    return list(set([e.scan for e in self]))

  def crystals(self):
    return list(set([e.crystal for e in self]))

  def imagesets(self):
    return list(set([e.imageset for e in self]))


class ExperimentManagerFactory(object):

  @staticmethod
  def from_datablock(datablock):
    pass

  @staticmethod
  def from_imageset_list(imageset_list):
    pass

  @staticmethod
  def from_imageset(imageset):
    pass





from dials.model.experiment import Beam, Detector, Goniometer, Scan, Crystal

if __name__ == '__main__':

  crystal0 = Crystal((1, 0, 0), (0, 1, 0), (0, 0, 1), 1)
  crystal1 = Crystal((1, 0, 0), (0, 1, 0), (0, 0, 1), 1)
  crystal2 = Crystal((1, 0, 0), (0, 1, 0), (0, 0, 1), 1)

  detector0 = Detector()
  detector1 = Detector()
  detector2 = detector0

  beam0 = Beam()
  beam1 = beam0
  beam2 = beam1

  expr0 = Experiment(beam=beam0, detector=detector0, crystal=crystal0)
  expr1 = Experiment(beam=beam1, detector=detector1, crystal=crystal1)
  expr2 = Experiment(beam=beam2, detector=detector2, crystal=crystal2)

  em = ExperimentManager()
  em.append(expr0)
  em.append(expr1)
  em.append(expr2)

  print "Unique Crystals"
  print em.crystals()

  print "Unique Beams"
  print em.beams()

  print "Unique Detectors"
  print em.detectors()

  print "Experiments with crystal:"
  print em.indices(crystal0)
  print em.indices(crystal1)
  print em.indices(crystal2)

  print "Experiments with beam:"
  print em.indices(beam0)
  print em.indices(beam1)
  print em.indices(beam2)

  print "Experiments with detector:"
  print em.indices(detector0)
  print em.indices(detector1)
  print em.indices(detector2)

  print "Experiments with detector/crystal combinations"
  print em.indices((detector0, crystal0))
  print em.indices((detector0, crystal1))
  print em.indices((detector0, crystal2))
  print em.indices((detector1, crystal0))
  print em.indices((detector1, crystal1))
  print em.indices((detector1, crystal2))
  print em.indices((detector2, crystal0))
  print em.indices((detector2, crystal1))
  print em.indices((detector2, crystal2))

  print "Replace model"
  em.replace(detector1, detector0)
  print em.detectors()

  print detector0 in em
  print detector1 in em

  em.extend(em)
  print em
