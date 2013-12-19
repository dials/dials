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

  def is_consistent(self):
    ''' If a scan is present, check that it makes sense with the imageset. '''
    from dxtbx.imageset2 import ImageSweep
    if self.scan:
      if isinstance(self.imageset, ImageSweep):
        if len(self.imageset) != self.scan.get_num_images():
          return False
        if self.imageset.get_array_range() != self.scan.get_array_range():
          return False
      elif self.imageset is not None:
        if (self.scan.get_num_images() != 1 or
            self.scan.get_oscillation()[1] != 0.0):
          return False
        if len(self.imageset.indices()) != 1:
          return False
        if self.imageset.indices()[0] != self.scan.get_array_range()[0]:
          return False
    return True

class ExperimentList(object):

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
      return ExperimentList(self._data[index])
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
    if isinstance(other, ExperimentList):
      self._data.extend(other._data)
    else:
      raise TypeError('expected type ExperimentList, got %s' % type(item))

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
    from collections import OrderedDict
    return OrderedDict([(e.beam, None) for e in self]).keys()

  def detectors(self):
    from collections import OrderedDict
    return OrderedDict([(e.detector, None) for e in self]).keys()

  def goniometers(self):
    from collections import OrderedDict
    return OrderedDict([(e.goniometer, None) for e in self]).keys()

  def scans(self):
    from collections import OrderedDict
    return OrderedDict([(e.scan, None) for e in self]).keys()

  def crystals(self):
    from collections import OrderedDict
    return OrderedDict([(e.crystal, None) for e in self]).keys()

  def imagesets(self):
    from collections import OrderedDict
    temp = OrderedDict([(e.imageset.reader(), i) for i, e in enumerate(self)
      if e is not None])
    return OrderedDict([(self[i].imageset.complete_set(), None)
      for i in temp.itervalues()]).keys()

  def is_consistent(self):
    return all([e.is_consistent() for e in self])

  def to_dict(self):
    from collections import OrderedDict
    from dxtbx.imageset2 import ImageSet, ImageSweep

    assert(self.is_consistent())

    # Get the list of unique models
    blist = self.beams()
    dlist = self.detectors()
    glist = self.goniometers()
    slist = self.scans()
    clist = self.crystals()
    ilist = self.imagesets()

    # Create the output dictionary
    result = OrderedDict()
    result['__id__'] = 'ExperimentList'
    result['experiment'] = []

    # Add the experiments to the dictionary
    for e in self:
      obj = OrderedDict()
      obj['__id__'] = 'Experiment'
      if e.beam:       obj['beam']       = blist.index(e.beam)
      if e.detector:   obj['detector']   = dlist.index(e.detector)
      if e.goniometer: obj['goniometer'] = glist.index(e.goniometer)
      if e.scan:       obj['scan']       = slist.index(e.scan)
      if e.crystal:    obj['crystal']    = clist.index(e.crystal)
      if e.imageset:
        obj['imageset'] = ilist.index(e.imageset)
        if e.scan is None and not isinstance(e.imageset, ImageSweep):
          if len(e.imageset) != len(e.imageset.complete_set()):
            obj['imageset'] = (obj['imageset'], e.imageset.indices())
      result['experiment'].append(obj)

    # Serialize all the imagesets
    result['imageset'] = []
    for imset in ilist:
      if isinstance(imset, ImageSweep):
        result['imageset'].append(OrderedDict([
          ('__id__', 'ImageSweep'),
          ('template', imset.get_template())]))
      elif isinstance(imset, ImageSet):
        result['imageset'].append(OrderedDict([
          ('__id__', 'ImageSet'),
          ('images', imset.paths())]))
      else:
        raise TypeError('expected ImageSet or ImageSweep, got %s' % type(imset))

    # Extract all the model dictionaries
    result['beam']       = [b.to_dict() for b in blist if b is not None]
    result['detector']   = [d.to_dict() for d in dlist if d is not None]
    result['goniometer'] = [g.to_dict() for g in glist if g is not None]
    result['scan']       = [s.to_dict() for s in slist if s is not None]
    result['crystal']    = [c.to_dict() for c in clist if c is not None]

    # Return the dictionary
    return result


class ExperimentListDict(object):

  def __init__(self, obj):
    from copy import deepcopy
    self._obj = deepcopy(obj)

  def decode(self):

    # Extract lists of models referenced by experiments
    self._blist = self._extract_models('beam')
    self._dlist = self._extract_models('detector')
    self._glist = self._extract_models('goniometer')
    self._slist = self._extract_models('scan')
    self._clist = self._extract_models('crystal')

    # Go through all the imagesets and make sure the dictionary
    # references by an index rather than a file path. Experiments
    # referencing the same imageset will get different objects
    # due to the fact that we can have different models
    self._ilist = self._extract_imagesets()

    # Extract all the experiments
    return self._extract_experiments()

  def _extract_models(self, name):

    # The from dict function
    from_dict = getattr(self, '_%s_from_dict' % name)

    # Extract all the model list
    mlist = self._obj.get(name, [])

    # Convert the model from dictionary to concreate
    # python class for the model.
    mlist = [from_dict(d) for d in mlist]

    # Dictionaries for file mappings
    mmap = {}

    # For each experiment, check the model is not specified by
    # a path, if it is then get the dictionary of the model
    # and insert it into the list. Replace the path reference
    # with an index
    for eobj in self._obj['experiment']:
      value = eobj.get(name, None)
      if value is None:
        continue
      elif isinstance(value, str):
        if value not in mmap:
          mmap[value] = len(mlist)
          mlist.append(from_dict(ExperimentListDict._from_file(value)))
        eobj[name] = mmap[value]
      elif not isinstance(value, int):
        raise TypeError('expected int or str, got %s' % type(value))

    # Return the model list
    return mlist

  def _extract_imagesets(self):

    # Extract all the model list
    mlist = self._obj.get('imageset', [])

    # Dictionaries for file mappings
    mmap = {}

    # For each experiment, check the imageset is not specified by
    # a path, if it is then get the dictionary of the imageset
    # and insert it into the list. Replace the path reference
    # with an index
    for eobj in self._obj['experiment']:
      value = eobj.get('imageset', None)
      if value is None:
        continue
      elif isinstance(value, str):
        if value not in mmap:
          mmap[value] = len(mlist)
          mlist.append(ExperimentListDict._from_file(value))
        eobj['imageset'] = mmap[value]
      elif not isinstance(value, int):
        raise TypeError('expected int or str, got %s' % type(value))

    # Return the model list
    return mlist

  def _extract_experiments(self):
    from dials.model.experiment.manager import ExperimentList

    # For every experiment, use the given input to create
    # a sensible experiment.
    el = ExperimentList()
    for eobj in self._obj['experiment']:
      el.append(self._create_experiment(
        ExperimentListDict.model_or_none(self._ilist, eobj, 'imageset'),
        ExperimentListDict.model_or_none(self._blist, eobj, 'beam'),
        ExperimentListDict.model_or_none(self._dlist, eobj, 'detector'),
        ExperimentListDict.model_or_none(self._glist, eobj, 'goniometer'),
        ExperimentListDict.model_or_none(self._slist, eobj, 'scan'),
        ExperimentListDict.model_or_none(self._clist, eobj, 'crystal')))

    # Return the experiment list
    return el

  def _create_experiment(self, imageset, beam, detector, goniometer, scan,
                         crystal):

    # Create the imageset from the input data
    if imageset is None:
      imageset = self._make_null()
    elif imageset['__id__'] == 'ImageSet':
      imageset = self._make_stills(imageset)
    elif imageset['__id__'] == 'ImageSweep':
      imageset = self._make_sweep(imageset, scan)

    # Fill in any models if they aren't already there
    if beam is None:
      beam = imageset.get_beam()
    if detector is None:
      detector = imageset.get_detector()
    if goniometer is None:
      goniometer = imageset.get_goniometer()
    if scan is None:
      scan = imageset.get_scan()

    # Return the experiment instance
    return Experiment(
      imageset=imageset,
      beam=beam,
      detector=detector,
      goniometer=goniometer,
      scan=scan,
      crystal=crystal
    )

  def _make_null(self):
    raise RuntimeError('NullSet not yet supported')

  def _make_stills(self, imageset):
    from dxtbx.imageset2 import ImageSetFactory
    return ImageSetFactory.make_imageset(imageset['images'])

  def _make_sweep(self, imageset, scan):
    from os.path import abspath, expanduser, expandvars
    from dxtbx.sweep_filenames import template_image_range
    from dxtbx.format.Registry import Registry
    from dxtbx.imageset2 import ImageSetFactory

    # Get the template format
    template = abspath(expanduser(expandvars(imageset['template'])))
    pfx = template.split('#')[0]
    sfx = template.split('#')[-1]
    template_format = '%s%%0%dd%s' % (pfx, template.count('#'), sfx)

    # Get the number of images (if no scan is given we'll try
    # to find all the images matching the template
    if scan is None:
      i0, i1 = template_image_range(template)
    else:
      i0, i1 = scan.get_image_range()

    # Get the format class from the first image
    format_class = Registry.find(template_format % i0)

    # Make a sweep from the input data
    return ImageSetFactory.make_sweep(template,
      list(range(i0, i1+1)), format_class)

  @staticmethod
  def model_or_none(mlist, eobj, name):
    index = eobj.get(name, None)
    if index is not None:
      return mlist[index]
    return None

  @staticmethod
  def _beam_from_dict(obj):
    from dxtbx.model import Beam
    return Beam.from_dict(obj)

  @staticmethod
  def _detector_from_dict(obj):
    from dxtbx.model import Detector, HierarchicalDetector
    if 'hierarchy' in obj:
      return HierarchicalDetector.from_dict(obj)
    else:
      return Detector.from_dict(obj)

  @staticmethod
  def _goniometer_from_dict(obj):
    from dxtbx.model import Goniometer
    return Goniometer.from_dict(obj)

  @staticmethod
  def _scan_from_dict(obj):
    from dxtbx.model import Scan
    return Scan.from_dict(obj)

  @staticmethod
  def _crystal_from_dict(obj):
    from dials.model.serialize import crystal
    return crystal.crystal_from_dict(obj)

  @staticmethod
  def _from_file(filename):
    from dxtbx.serialize.load import _decode_dict
    from os.path import expanduser, expandvars, abspath
    import json
    filename = abspath(expanduser(expandvars(filename)))
    try:
      with open(filename, 'r') as infile:
        return json.loads(infile.read(), object_hook=_decode_dict)
    except IOError, e:
      raise IOError('unable to read file, %s' % filename)


class ExperimentListFactory(object):

  @staticmethod
  def from_args(args, verbose=False, unhandled=None):
    ''' Try to load datablocks from any recognised format. '''
    from dxtbx.datablock import DataBlockFactory

    # Create a list for unhandled arguments
    if unhandled is None:
      unhandled = []
    unhandled1 = []
    unhandled2 = []

    # First try as image files
    experiments = ExperimentListFactory.from_datablock(
      DataBlockFactory.from_args(args, verbose, unhandled1))

    # Now try as JSON files
    if len(unhandled1) > 0:
      for filename in unhandled1:
        try:
          experiments.extend(ExperimentListFactory.from_json_file(filename))
          if verbose: print 'Loaded experiment(s) from %s' % filename
        except Exception, e:
          unhandled2.append(filename)

    # Now try as pickle files
    if len(unhandled2) > 0:
      for filename in unhandled2:
        try:
          experiments.extend(ExperimentListFactory.from_pickle_file(filename))
          if verbose: print 'Loaded experiments(s) from %s' % filename
        except Exception:
          unhandled.append(filename)

    # Return the experiments
    return experiments

  @staticmethod
  def from_datablock(datablock):

    # Initialise the experiment list
    experiments = ExperimentList()

    # If we have a list, loop through
    if isinstance(datablock, list):
      for db in datablock:
        experiments.extend(ExperimentListFactory.from_datablock(db))
      return experiments

    # Extract the sweeps
    for sweep in datablock.extract_sweeps():
      experiments.append(Experiment(
        imageset=sweep,
        beam=sweep.get_beam(),
        detector=sweep.get_detector(),
        goniometer=sweep.get_goniometer(),
        scan=sweep.get_scan()))

    # Extract the stills
    stills = datablock.extract_stills()
    for i in range(len(stills)):
      still = stills[i:i+1]
      experiments.append(Experiment(
        imageset=still,
        beam=still.get_beam(),
        detector=still.get_detector(),
        goniometer=still.get_goniometer(),
        scan=still.get_scan()))

    # Return the experiments
    return experiments

  @staticmethod
  def from_dict(obj):
    return ExperimentListDict(obj).decode()

  @staticmethod
  def from_json(text):
    from dxtbx.serialize.load import _decode_dict
    import json
    return ExperimentListFactory.from_dict(
      json.loads(text, object_hook=_decode_dict))

  @staticmethod
  def from_json_file(filename):
    with open(filename, 'r') as infile:
      return ExperimentListFactory.from_json(infile.read())

  @staticmethod
  def from_pickle_file(filename):
    ''' Decode an experiment list from a pickle file. '''
    import cPickle as pickle
    with open(filename, 'rb') as infile:
      obj = pickle.load(infile)
      assert(isinstance(obj, ExperimentList))
      return obj


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

  em = ExperimentList()
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
