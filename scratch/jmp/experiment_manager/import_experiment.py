from dials.model.experiment.manager import Experiment, ExperimentManager

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
    from dials.model.experiment.manager import ExperimentManager

    # For every experiment, use the given input to create
    # a sensible experiment.
    el = ExperimentManager()
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
    from dxtbx.imageset import ImageSetFactory
    return ImageSetFactory.make_imageset(imageset['images'])

  def _make_sweep(self, imageset, scan):
    from os.path import abspath, expanduser, expandvars
    from dxtbx.sweep_filenames import template_image_range
    from dxtbx.format.Registry import Registry
    from dxtbx.imageset import ImageSetFactory

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
    from cctbx.crystal.crystal_model.serialize import crystal_from_dict
    return crystal_from_dict(obj)

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


def from_json(text):
  from dxtbx.serialize.load import _decode_dict
  import json
  return ExperimentListDict(json.loads(
    text, object_hook=_decode_dict)).decode()

def from_json_file(filename):
  with open(filename, 'r') as infile:
    return from_json(infile.read())


if __name__ == '__main__':
  import os
  import sys
  import libtbx.load_env
  try:
    dials_regression = libtbx.env.dist_path('dials_regression')
  except KeyError, e:
    print 'FAIL: dials_regression not configured'
    exit(0)
  os.environ['DIALS_REGRESSION'] = dials_regression


  num = 1
  if len(sys.argv) == 2:
    num = sys.argv[1]

  filename = os.path.join(dials_regression, 'experiment_test_data',
                          'experiment_{0}.json'.format(num))


  el = from_json_file(filename)
  for e in el:
    print e.imageset
    print e.beam
    print e.detector
    print e.goniometer
    print e.scan
    print e.crystal
  #print len(el)
  #print el[0].imageset
  #print el[0].beam
  #print el[0].detector
  #print el[0].goniometer
  #print el[0].scan
  #print el[0].crystal
