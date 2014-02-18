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
  def from_predictions(exlist, force_static=False):
    ''' Construct a reflection table from predictions. '''
    from dials.algorithms.spot_prediction.reflection_predictor \
      import ReflectionPredictor
    predict = ReflectionPredictor(exlist, force_static=force_static)
    return predict()

  @staticmethod
  def from_observations(datablocks, params):
    ''' Construct a reflection table from observations. '''
    from dials.algorithms.peak_finding.spotfinder_factory \
      import SpotFinderFactory

    # Ensure we have a data block
    if len(datablocks) != 1:
      raise RuntimeError('only 1 datablock can be processed at a time')

    # Get the integrator from the input parameters
    print 'Configuring spot finder from input parameters'
    find_spots = SpotFinderFactory.from_parameters(params)

    # Find the spots
    return find_spots(datablocks[0])

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

  def compute_bbox(self, experiment, nsigma, sigma_d=None, sigma_m=None):
    ''' Compute the bounding boxes. '''
    from dials.algorithms.shoebox import BBoxCalculator
    from dials.util.command_line import Command

    # Get the beam divergence and mosaicity
    if sigma_d is None or sigma_m is None:
      sigma_d = experiment.beam.get_sigma_divergence(deg=False)
      sigma_m = experiment.crystal.get_mosaicity(deg=False)

    # Create the bbox calculator
    calculate = BBoxCalculator(
      experiment.beam, experiment.detector,
      experiment.goniometer, experiment.scan,
      nsigma * sigma_d,
      nsigma * sigma_m)

    # Calculate the bounding boxes of all the reflections
    Command.start('Calculating bounding boxes')
    self['bbox'] = calculate(
      self['s1'],
      self['xyzcal.mm'].parts()[2],
      self['panel'])
    Command.end('Calculated {0} bounding boxes'.format(len(self)))

  def compute_background(self, experiment, parameters):
    ''' Helper function to compute the background. '''
    from dials.algorithms.background.background_factory import BackgroundFactory
    function = BackgroundFactory.from_parameters(parameters)
    function(experiment, self)

  def compute_centroid(self, experiment, parameters):
    ''' Helper function to compute the centroid. '''
    from dials.algorithms.centroid.centroid_factory import CentroidFactory
    function = CentroidFactory.from_parameters(parameters)
    function(experiment, self)

  def compute_intensity(self, experiment, parameters):
    ''' Helper function to compute the intensity. '''
    from dials.algorithms.integration.integrator import IntensityFactory
    function = IntensityFactory.from_parameters(parameters)
    function(experiment, self)

  def correct_intensity(self, experiment):
    ''' Helper function to correct the intensity. '''
    from dials.algorithms.integration.lp_correction import correct_intensity
    correct_intensity(experiment, self)

  def integrate(self, experiment, parameters, save_profiles=False):
    ''' Helper function to integrate reflections. '''
    self.compute_background(experiment, parameters)
    self.compute_centroid(experiment, parameters)
    self.compute_intensity(experiment, parameters)
    self.correct_intensity(experiment)
    if save_profiles == False:
      del self['shoebox']
