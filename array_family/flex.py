from __future__ import division
import boost.python
from cctbx.array_family.flex import *
from dials.model import data
from dials_array_family_flex_ext import *

# Set the 'real' type to either float or double
if get_real_type() == "float":
  real = float
elif get_real_type() == "double":
  real = double
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
