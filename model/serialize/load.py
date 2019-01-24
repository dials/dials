from __future__ import absolute_import, division, print_function

# Import to give access from here
from dxtbx.serialize.load import imageset as sweep # implicit import

import six.moves.cPickle as pickle

def reflections(infile):
  '''
  Load the given reflection file.

  :params infile: The input filename or file object
  :returns: The reflection list
  '''

  # If the input is a string then open and read from that file
  if isinstance(infile, str):
    with open(infile, 'rb') as infile:
      return pickle.load(infile)

  # Otherwise assume the input is a file and read from it
  else:
    return pickle.load(infile)

def reference(infile):
  '''
  Load the given reference profile file.

  :params infile: The input filename or file object
  :returns: The reference list
  '''

  # If the input is a string then open and read from that file
  if isinstance(infile, str):
    with open(infile, 'rb') as infile:
      return pickle.load(infile)

  # Otherwise assume the input is a file and read from it
  else:
    return pickle.load(infile)
