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
  from dials.extensions import SimpleBackgroundExt
  return strategy(SimpleBackgroundExt)

def default_intensity_algorithm():
  ''' Get the default intensity algorithm. '''
  from dials.extensions import SummationIntegrationExt
  return strategy(SummationIntegrationExt)

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
  def from_predictions(experiment,
                       dmin=None,
                       dmax=None,
                       margin=1,
                       force_static=False):
    ''' Construct a reflection table from predictions. '''
    from dials.algorithms.spot_prediction.reflection_predictor \
      import ReflectionPredictor
    predict = ReflectionPredictor(
      experiment,
      dmin=dmin,
      dmax=dmax,
      margin=margin,
      force_static=force_static)
    return predict()

  @staticmethod
  def from_predictions_multi(experiments,
                             dmin=None,
                             dmax=None,
                             margin=1,
                             force_static=False):
    ''' Construct a reflection table from predictions. '''
    from scitbx.array_family import flex
    result = reflection_table()
    for i, e in enumerate(experiments):
      rlist = reflection_table.from_predictions(
        e,
        dmin=dmin,
        dmax=dmax,
        margin=margin,
        force_static=force_static)
      rlist['id'] = flex.size_t(len(rlist), i)
      result.extend(rlist)
    return result

  @staticmethod
  def from_observations(datablock, params=None):
    ''' Construct a reflection table from observations. '''
    from dials.algorithms.peak_finding.spotfinder_factory \
      import SpotFinderFactory
    from logging import info

    # Get the integrator from the input parameters
    info('Configuring spot finder from input parameters')
    find_spots = SpotFinderFactory.from_parameters(params)

    # Find the spots
    return find_spots(datablock)

  @staticmethod
  def from_pickle(filename):
    ''' Read the reflection table from pickle file. '''
    import cPickle as pickle
    with open(filename, 'rb') as infile:
      result = pickle.load(infile)
      assert(isinstance(result, reflection_table))
      return result

  @staticmethod
  def from_h5(filename):
    ''' Read the reflections table from a HDF5 file. '''
    from dials.util.nexus import NexusFile
    handle = NexusFile(filename, 'r')
    self = handle.get_reflections()
    handle.close()
    return self

  @staticmethod
  def empty_standard(nrows):
    ''' Create an empty table of specified number of rows with most of the
    standard keys'''

    assert nrows > 0
    table = reflection_table(nrows)

    # General properties
    table['flags'] = flex.size_t(nrows, 0)
    table['id'] = flex.double(nrows, 0)
    table['panel'] = flex.size_t(nrows, 0)

    # Predicted properties
    table['miller_index'] = flex.miller_index(nrows)
    table['entering'] = flex.bool(nrows)
    table['s1'] = flex.vec3_double(nrows, (0, 0, 0))
    table['xyzcal.mm'] = flex.vec3_double(nrows, (0, 0, 0))
    table['xyzcal.px'] = flex.vec3_double(nrows, (0, 0, 0))
    #table['ub_matrix'] = flex.mat3_double(nrows, (0, 0, 0, 0, 0, 0, 0, 0, 0))

    # Observed properties
    table['xyzobs.px.value'] = flex.vec3_double(nrows, (0, 0, 0))
    table['xyzobs.px.variance'] = flex.vec3_double(nrows, (0, 0, 0))
    table['xyzobs.mm.value'] = flex.vec3_double(nrows, (0, 0, 0))
    table['xyzobs.mm.variance'] = flex.vec3_double(nrows, (0, 0, 0))
    table['rlp'] = flex.vec3_double(nrows, (0, 0, 0))
    table['intensity.sum.value'] = flex.double(nrows, 0)
    table['intensity.sum.variance'] = flex.double(nrows, 0)
    table['intensity.prf.value'] = flex.double(nrows, 0)
    table['intensity.prf.variance'] = flex.double(nrows, 0)
    table['lp'] = flex.double(nrows, 0)
    table['profile.correlation'] = flex.double(nrows, 0)

    return table

  @staticmethod
  def plot(table, detector, key):
    ''' Plot a reflection table using matplotlib '''
    from matplotlib import pyplot as plt
    from matplotlib.patches import Polygon
    fig=plt.figure()
    ax=fig.add_subplot(111,aspect='equal')
    spots = table[key]
    if 'px' in key:
      spots = [detector[table['panel'][i]].get_pixel_lab_coord(spots[i][0:2]) for i in xrange(len(spots))]
    else:
      assert 'mm' in key
      spots = [detector[table['panel'][i]].get_lab_coord(spots[i][0:2]) for i in xrange(len(spots))]

    min_f = max_f = min_s = max_s = 0

    for i, panel in enumerate(detector):
      fs, ss = panel.get_image_size()
      p0 = panel.get_pixel_lab_coord((0,0))
      p1 = panel.get_pixel_lab_coord((fs-1,0))
      p2 = panel.get_pixel_lab_coord((fs-1,ss-1))
      p3 = panel.get_pixel_lab_coord((0,ss-1))
      p = Polygon((p0[0:2],p1[0:2],p2[0:2],p3[0:2]), closed=True, color='green', fill=False, hatch='/')

      if p.xy[:,0].min() < min_f: min_f = p.xy[:,0].min()
      if p.xy[:,0].max() > max_f: max_f = p.xy[:,0].max()
      if p.xy[:,1].min() < min_s: min_s = p.xy[:,1].min()
      if p.xy[:,1].max() > max_s: max_s = p.xy[:,1].max()

      ax.add_patch(p)

    ax.set_xlim((min_f-10,max_f+10))
    ax.set_ylim((min_s-10,max_s+10))
    plt.scatter([s[0] for s in spots],
                [s[1] for s in spots],c='blue',linewidth=0)
    plt.show()


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

  def match(self, other):
    ''' Match reflections with another set of reflections. '''
    from dials.algorithms.peak_finding.spot_matcher import SpotMatcher
    match = SpotMatcher(max_separation=1)
    oind, sind = match(other, self)
    return sind, oind

  def match_with_reference(self, other):
    ''' Match reflections with another set of reflections. '''
    from logging import info
    info("Matching reference spots with predicted reflections")
    sind, oind = self.match(other)
    h1 = self.select(sind)['miller_index']
    h2 = other.select(oind)['miller_index']
    mask = (h1 == h2)
    self.set_flags(sind.select(mask), self.flags.reference_spot)
    info("Matched %d reference spots with predicted reflections" %
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
    from dials.algorithms.profile_model.gaussian_rs import zeta_factor
    m2 = experiment.goniometer.get_rotation_axis()
    s0 = experiment.beam.get_s0()
    self['zeta'] = zeta_factor(m2, s0, self['s1'])
    return self['zeta']

  def compute_zeta_multi(self, experiments):
    ''' Compute zeta for each reflection. '''
    from dials.algorithms.profile_model.gaussian_rs import zeta_factor
    m2 = flex.vec3_double(len(experiments))
    s0 = flex.vec3_double(len(experiments))
    for i, e in enumerate(experiments):
      m2[i] = e.goniometer.get_rotation_axis()
      s0[i] = e.beam.get_s0()
    self['zeta'] = zeta_factor(m2, s0, self['s1'], self['id'])
    return self['zeta']

  def compute_d_single(self, experiment):
    ''' Compute the resolution for each reflection. '''
    from dials.array_family import flex
    uc = flex.unit_cell(1)
    uc[0] = experiment.crystal.get_unit_cell()
    self['d'] = uc.d(self['miller_index'], flex.size_t(len(self), 0))
    return self['d']

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
    return profile_model.compute_bbox(
      experiments,
      self,
      sigma_b_multiplier=sigma_b_multiplier)

  def compute_partiality(self, experiments, profile_model):
    ''' Compute the reflection partiality. '''
    return profile_model.compute_partiality(experiments, self)

  def compute_background(self, experiments):
    ''' Helper function to compute the background. '''
    self._background_algorithm(experiments).compute_background(self)

  def compute_centroid(self, experiments):
    ''' Helper function to compute the centroid. '''
    self._centroid_algorithm(experiments).compute_centroid(self)

  def compute_intensity(self, experiments, profile_model):
    ''' Helper function to compute the intensity. '''
    self._intensity_algorithm(experiments, profile_model).compute_intensity(self)

  def compute_summed_intensity(self):
    ''' Compute intensity via summation integration. '''
    from dials.algorithms.integration.sum import IntegrationAlgorithm
    algorithm = IntegrationAlgorithm()
    algorithm(self)

  def compute_corrections(self, experiments):
    ''' Helper function to correct the intensity. '''
    from dials.algorithms.integration import Corrections, CorrectionsMulti
    from logging import info
    info("Calculating lp correction")
    compute = CorrectionsMulti()
    for experiment in experiments:
      compute.append(Corrections(
        experiment.beam,
        experiment.goniometer))
    lp = compute.lp(self['id'], self['s1'])
    self['lp'] = lp
    info("Calculated lp correction")
    return lp

  def integrate(self, experiments, profile_model):
    ''' Helper function to integrate reflections. '''
    self.compute_background(experiments)
    self.compute_centroid(experiments)
    self.compute_summed_intensity()
    self.compute_intensity(experiments, profile_model)

  def compute_mask(self, experiments, profile_model):
    ''' Apply a mask to the shoeboxes. '''
    from dials.algorithms.shoebox import Masker3DProfile
    assert(len(experiments) == 1)
    mask_profiles = Masker3DProfile(experiments, profile_model)
    mask_profiles(self, None)

  def extract_shoeboxes(self, imageset, mask=None):
    ''' Helper function to read a load of shoebox data. '''
    from dials.model.data import make_image
    from logging import info
    from time import time
    assert("shoebox" in self)
    detector = imageset.get_detector()
    try:
      frame0, frame1 = imageset.get_array_range()
    except Exception:
      frame0, frame1 = (0, len(imageset))
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
    info("Reading images")
    read_time = 0
    extract_time = 0
    for i in range(len(imageset)):
      st = time()
      image = imageset[i]
      read_time += time() - st
      if not isinstance(image, tuple):
        image = (image,)
      st = time()
      extractor.next(make_image(image, mask))
      extract_time += time() - st
      del image
    assert(extractor.finished())
    return read_time, extract_time

  def is_overloaded(self, experiments):
    ''' Check if the shoebox contains overloaded pixels. '''
    from dials.algorithms.shoebox import OverloadChecker
    assert('shoebox' in self)
    assert('id' in self)
    checker = OverloadChecker()
    for experiment in experiments:
      checker.add(flex.double(
        (p.get_trusted_range()[1]
          for p in experiment.detector)))
    result = checker(self['id'], self['shoebox'])
    self.set_flags(result, self.flags.overloaded)
    return result
