from __future__ import absolute_import, division, print_function

# Import to give access from here
from dxtbx.serialize.dump import imageset as sweep # implicit import
from dxtbx.serialize.dump import imageset_to_string as sweep_to_string # implicit import
from dxtbx.serialize.dump import datablock # implicit import

import six.moves.cPickle as pickle

def reflections(obj, outfile):
  '''
  Dump the given object to file

  :param obj: The reflection list to dump
  :param outfile: The output file name or file object
  '''
  if isinstance(outfile, str):
    with open(outfile, 'wb') as outfile:
      pickle.dump(obj, outfile, pickle.HIGHEST_PROTOCOL)

  # Otherwise assume the input is a file and write to it
  else:
    pickle.dump(obj, outfile, pickle.HIGHEST_PROTOCOL)

def reference(obj, outfile):
  '''
  Dump the given object to file

  :param obj: The reference list to dump
  :param outfile: The output file name or file object

  '''

  if isinstance(outfile, str):
    with open(outfile, 'wb') as outfile:
      pickle.dump(obj, outfile, pickle.HIGHEST_PROTOCOL)

  # Otherwise assume the input is a file and write to it
  else:
    pickle.dump(obj, outfile, pickle.HIGHEST_PROTOCOL)
