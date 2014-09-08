#
# flex.py
#
#  Copyright (C) 2013 Diamond Light Source
#
#  Author: James Parkhurst
#
#  This code is distributed under the BSD license, a copy of which is
#  included in the root directory of this package.
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

def strategy(cls, params=None):
  ''' Wrap a class that takes params and experiments as a strategy. '''
  from functools import wraps
  @wraps(cls)
  def call(self, *args):
    return cls(params, *args)
  return call

def default_background_algorithm():
  ''' Get the default background algorithm. '''
  from dials.extensions import GeneralBackgroundExt
  return strategy(GeneralBackgroundExt)

def default_intensity_algorithm():
  ''' Get the default intensity algorithm. '''
  from dials.extensions import Summation3dIntegrationExt
  return strategy(Summation3dIntegrationExt)

def default_centroid_algorithm():
  ''' Get the default centroid algorithm. '''
  from dials.extensions import SimpleCentroidExt
  return strategy(SimpleCentroidExt)


class reflection_table_aux(boost.python.injector, reflection_table):
  ''' An injector class to add additional methods to the reflection table. '''

  # Set the default algorithms. These are set as class variables so that if they
  # are changed in the class, all new instances of reflection table will have
  # the modified algorithms. If these are modified on the instance level, then
  # only the instance will have the modified algorithms and new instances will
  # have the defaults
  _background_algorithm = default_background_algorithm()
  _intensity_algorithm = default_intensity_algorithm()
  _centroid_algorithm = default_centroid_algorithm()

  @staticmethod
  def from_predictions(experiment, force_static=False, dmin=None):
    ''' Construct a reflection table from predictions. '''
    from dials.algorithms.spot_prediction.reflection_predictor \
      import ReflectionPredictor
    predict = ReflectionPredictor(
      experiment, force_static=force_static, dmin=dmin)
    return predict()

  @staticmethod
  def from_predictions_multi(experiments, force_static=False, dmin=None):
    ''' Construct a reflection table from predictions. '''
    from scitbx.array_family import flex
    result = reflection_table()
    for i, e in enumerate(experiments):
      rlist = reflection_table.from_predictions(e, force_static, dmin)
      rlist['id'] = flex.size_t(len(rlist), i)
      result.extend(rlist)
    return result

  @staticmethod
  def from_observations(datablock, params=None):
    ''' Construct a reflection table from observations. '''
    from dials.algorithms.peak_finding.spotfinder_factory \
      import SpotFinderFactory

    # Get the integrator from the input parameters
    print 'Configuring spot finder from input parameters'
    find_spots = SpotFinderFactory.from_parameters(params)

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

  def copy(self):
    ''' Copy everything. '''
    from scitbx.array_family import flex
    return self.select(flex.bool(len(self), True))

  def sort(self, name, reverse=False):
    ''' Sort the reflection table by a key. '''
    import __builtin__
    column = self[name]
    indices = __builtin__.sorted(
      range(len(self)),
      key=lambda x: column[x],
      reverse=reverse)
    self.reorder(flex.size_t(indices))

  def split_by_experiment_id(self):
    ''' Split the reflection table into multiple tables by experiment id. '''
    from scitbx.array_family import flex
    temp = self.select(flex.bool(len(self), True))
    result = []
    i = 0
    while (True):
      mask = temp['id'] == i
      new_list = temp.select(mask)
      temp.del_selected(mask)
      if len(new_list) > 0:
        result.append(new_list)
      if len(temp) == 0:
        break
      i += 1
    return result

  def match(self, other):
    ''' Match reflections with another set of reflections. '''
    from dials.algorithms.peak_finding.spot_matcher import SpotMatcher
    match = SpotMatcher(max_separation=1)
    oind, sind = match(other, self)
    return sind, oind

  def match_with_reference(self, other):
    ''' Match reflections with another set of reflections. '''
    from dials.util.command_line import Command
    Command.start("Matching reference spots with predicted reflections")
    sind, oind = self.match(other)
    h1 = self.select(sind)['miller_index']
    h2 = other.select(oind)['miller_index']
    mask = (h1 == h2)
    self.set_flags(sind.select(mask), self.flags.reference_spot)
    Command.end("Matched %d reference spots with predicted reflections" %
                mask.count(True))

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

  def compute_bbox(self, experiments, profile_model, sigma_b_multiplier=2.0):
    ''' Compute the bounding boxes. '''
    from dials.util.command_line import Command
    profile_model.compute_bbox(experiments, self, sigma_b_multiplier)

  def compute_partiality(self, experiments, profile_model):
    ''' Compute the reflection partiality. '''
    profile_model.compute_partiality(experiments, self)

  def compute_background(self, experiments):
    ''' Helper function to compute the background. '''
    self._background_algorithm(experiments).compute_background(self)

  def compute_centroid(self, experiments):
    ''' Helper function to compute the centroid. '''
    self._centroid_algorithm(experiments).compute_centroid(self)

  def compute_intensity(self, experiments, profile_model):
    ''' Helper function to compute the intensity. '''
    self._intensity_algorithm(experiments, profile_model).compute_intensity(self)

  def compute_corrections(self, experiments):
    ''' Helper function to correct the intensity. '''
    from dials.algorithms.integration import Corrections, CorrectionsMulti
    from dials.util.command_line import Command
    Command.start("Calculating lp correction")
    compute = CorrectionsMulti()
    for experiment in experiments:
      compute.append(Corrections(
        experiment.beam,
        experiment.goniometer))
    lp = compute.lp(self['id'], self['s1'])
    self['lp'] = lp
    Command.end("Calculated lp correction")
    return lp

  def integrate(self, experiments, profile_model):
    ''' Helper function to integrate reflections. '''
    self.compute_background(experiments)
    self.compute_centroid(experiments)
    self.compute_intensity(experiments, profile_model)

  def compute_mask(self, experiments, profile_model):
    ''' Apply a mask to the shoeboxes. '''
    from dials.algorithms.shoebox import Masker3DProfile
    from math import pi
    assert(len(experiments) == 1)
    mask_profiles = Masker3DProfile(experiments, profile_model)
    mask_profiles(self, None)

  def extract_shoeboxes(self, imageset, mask=None):
    ''' Helper function to read a load of shoebox data. '''
    from dials.model.data import Image
    import sys
    from time import time
    assert("shoebox" in self)
    detector = imageset.get_detector()
    frame0, frame1 = imageset.get_array_range()
    extractor = ShoeboxExtractor(self, len(detector), frame0, frame1)
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
