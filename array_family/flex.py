from __future__ import division
import boost.python
from dials.model import data
from dials_array_family_flex_ext import *
from cctbx.array_family.flex import *
from cctbx.array_family import flex

# Set the 'real' type to either float or double
if get_real_type() == "float":
  real = flex.float
elif get_real_type() == "double":
  real = flex.double
else:
  raise TypeError('unknown "real" type')


class reflection_table_aux(boost.python.injector, reflection_table):
  ''' An injector class to add additional methods to the reflection table. '''

  @staticmethod
  def from_predictions(experiment, force_static=False, dmin=None):
    ''' Construct a reflection table from predictions. '''
    from dials.algorithms.spot_prediction.reflection_predictor \
      import ReflectionPredictor
    predict = ReflectionPredictor(
      experiment, force_static=force_static, dmin=dmin)
    return predict()

  @staticmethod
  def from_observations(datablock):
    ''' Construct a reflection table from observations. '''
    from dials.algorithms.peak_finding.spotfinder_factory \
      import SpotFinderFactory

    # Get the integrator from the input parameters
    print 'Configuring spot finder from input parameters'
    from dials.framework.registry import Registry
    find_spots = SpotFinderFactory.from_parameters(Registry().params())

    # Find the spots
    return find_spots(datablock)

  @staticmethod
  def from_pickle(filename):
    ''' Read the reflection table from pickle file. '''
    import cPickle as pickle
    with open(filename, 'rb') as infile:
      return pickle.load(infile)

  @staticmethod
  def from_h5(filename):
    ''' Read the reflections table from a HDF5 file. '''
    from dials.util.nexus import NexusFile
    handle = NexusFile(filename, 'r')
    self = handle.get_reflections()
    handle.close()
    return self

  def as_pickle(self, filename):
    ''' Write the reflection table as a pickle file. '''
    import cPickle as pickle
    with open(filename, 'wb') as outfile:
      pickle.dump(self, outfile, protocol=pickle.HIGHEST_PROTOCOL)

  def as_h5(self, filename):
    ''' Write the reflection table as a HDF5 file. '''
    from dials.util.nexus import NexusFile
    handle = NexusFile(filename, 'w')
    handle.set_reflections(self)
    handle.close()

  def sort(self, name, reverse=False):
    ''' Sort the reflection table by a key. '''
    import __builtin__
    column = self[name]
    indices = __builtin__.sorted(
      range(len(self)),
      key=lambda x: column[x],
      reverse=reverse)
    self.reorder(flex.size_t(indices))

  #def is_bbox_inside_image_range(self, experiment):
    #''' Check if bbox is within image range. '''
    #from dials.algorithms import filtering
    #assert(len(experiment.detector) == 1)
    #return filtering.is_bbox_outside_image_range(
      #self['bbox'],
      #experiment.detector[0].get_image_size()[::-1],
      #experiment.scan.get_array_range()) != True

  def compute_zeta(self, experiment):
    ''' Compute zeta for each reflection. '''
    from dials.algorithms.reflection_basis import zeta_factor
    m2 = experiment.goniometer.get_rotation_axis()
    s0 = experiment.beam.get_s0()
    self['zeta'] = zeta_factor(m2, s0, self['s1'])
    return self['zeta']

  def compute_zeta_multi(self, experiments):
    ''' Compute zeta for each reflection. '''
    from dials.algorithms.reflection_basis import zeta_factor
    m2 = flex.vec3_double(len(experiments))
    s0 = flex.vec3_double(len(experiments))
    for i, e in enumerate(experiments):
      m2[i] = e.goniometer.get_rotation_axis()
      s0[i] = e.beam.get_s0()
    self['zeta'] = zeta_factor(m2, s0, self['s1'], self['id'])
    return self['zeta']

  def compute_d(self, experiments):
    ''' Compute the resolution for each reflection. '''
    from dials.array_family import flex
    uc = flex.unit_cell(len(experiments))
    for i, e in enumerate(experiments):
      uc[i] = e.crystal.get_unit_cell()
    self['d'] = uc.d(self['miller_index'], self['id'])
    return self['d']


  def compute_bbox(self, experiment, nsigma, sigma_d=None, sigma_m=None,
                   sigma_d_multiplier=2.0):
    ''' Compute the bounding boxes. '''
    from dials.algorithms.shoebox import BBoxCalculator
    from dials.util.command_line import Command
    from dials.framework.registry import Registry
    from math import pi

    # Get the beam divergence and mosaicity
    if sigma_d is None or sigma_m is None:
      registry = Registry()
      sigma_d = registry.params().integration.shoebox.sigma_b * pi / 180.0
      sigma_m = registry.params().integration.shoebox.sigma_m * pi / 180.0

    # Create the bbox calculator
    # sigma_d_multiplier so as to include a region of background pixels
    # in the shoebox
    Command.start('Calculating bounding boxes')
    calculate = BBoxCalculator(
      experiment.beam, experiment.detector,
      experiment.goniometer, experiment.scan,
      sigma_d_multiplier * nsigma * sigma_d,
      nsigma * sigma_m)

    # Calculate the bounding boxes of all the reflections
    self['bbox'] = calculate(
      self['s1'],
      self['xyzcal.mm'].parts()[2],
      self['panel'])
    Command.end('Calculated {0} bounding boxes'.format(len(self)))

  def compute_background(self, experiment):
    ''' Helper function to compute the background. '''
    from dials.framework.registry import init_ext
    init_ext("integration.background", experiment).compute_background(self)

  def compute_centroid(self, experiment):
    ''' Helper function to compute the centroid. '''
    from dials.framework.registry import init_ext
    init_ext("integration.centroid", experiment).compute_centroid(self)

  def compute_intensity(self, experiment):
    ''' Helper function to compute the intensity. '''
    from dials.framework.registry import init_ext
    init_ext("integration.intensity", experiment).compute_intensity(self)

  def correct_intensity(self, experiment):
    ''' Helper function to correct the intensity. '''
    if experiment.goniometer is not None:
      from dials.algorithms.integration.lp_correction import correct_intensity
      correct_intensity(experiment, self)

  def integrate(self, experiment, save_profiles=False):
    ''' Helper function to integrate reflections. '''
    self.compute_background(experiment)
    self.compute_centroid(experiment)
    self.compute_intensity(experiment)
    self.correct_intensity(experiment)

  def fill_shoeboxes(self, imageset, mask=None):
    ''' Helper function to read a load of shoebox data. '''
    from dials.model.serialize import SimpleShoeboxExtractor
    from dials.model.data import Image
    import sys
    from time import time
    assert("shoebox" in self)
    detector = imageset.get_detector()
    frame0, frame1 = imageset.get_array_range()
    extractor = SimpleShoeboxExtractor(
      self['shoebox'],
      frame0,
      frame1,
      len(detector))
    if mask is None:
      image = imageset[0]
      if not isinstance(image, tuple):
        image = (image,)
      mask = []
      for i in range(len(image)):
        tr = detector[i].get_trusted_range()
        mask.append(image[i].as_double() > tr[0])
      mask = tuple(mask)
    sys.stdout.write("Reading images: ")
    read_time = 0
    extract_time = 0
    for i in range(len(imageset)):
      st = time()
      image = imageset[i]
      read_time += time() - st
      if not isinstance(image, tuple):
        image = (image,)
      st = time()
      extractor.next(Image(image, mask))
      extract_time += time() - st
      sys.stdout.write(".")
      sys.stdout.flush()
      del image
    sys.stdout.write("\n")
    sys.stdout.flush()
    assert(extractor.finished())
    return read_time, extract_time

