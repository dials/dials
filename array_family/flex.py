from __future__ import division
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


@staticmethod
def reflection_table_from_predictions(exlist):
  ''' Construct a reflection table from predictions. '''
  from dials.algorithms.integration import ReflectionPredictor
  from dials.array_family import flex
  predict = ReflectionPredictor()
  result = flex.reflection_table()
  for idx, ex in enumerate(exlist):
    rlist = predict(ex.imageset, ex.crystal)
    rtable = rlist.to_table()
    rtable['id'] = flex.size_t(len(rlist), idx)
    result.extend(rtable)
  return result

def reflection_table_as_pickle(self, filename):
  ''' Write the reflection table as a pickle file. '''
  import cPickle as pickle
  with open(filename, 'wb') as outfile:
    pickle.dump(self, outfile, protocol=pickle.HIGHEST_PROTOCOL)

reflection_table.from_predictions = reflection_table_from_predictions
reflection_table.as_pickle = reflection_table_as_pickle
