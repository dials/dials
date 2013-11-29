#!/usr/bin/env python
#
# dials.model.serialize.load.py
#
#  Copyright (C) 2013 Diamond Light Source
#
#  Author: James Parkhurst
#
#  This code is distributed under the BSD license, a copy of which is
#  included in the root directory of this package.

from __future__ import division

# Import to give access from here
from dxtbx.serialize.load import imageset as sweep

# FIXME These are only temporary - they should be removed ASAP
from cctbx.crystal.crystal_model.serialize import crystal_from_string
from cctbx.crystal.crystal_model.serialize import load_crystal as crystal

def reflections(infile):
  ''' Load the given reflection file.

  Params:
      infile The input filename or file object

  Returns:
      The reflection list

  '''
  import cPickle as pickle

  # If the input is a string then open and read from that file
  if isinstance(infile, str):
    with open(infile, 'rb') as infile:
      return pickle.load(infile)

  # Otherwise assume the input is a file and read from it
  else:
    return pickle.load(infile)

def reference(infile):
  ''' Load the given reference profile file.

  Params:
      infile The input filename or file object

  Returns:
      The reflection list

  '''
  import cPickle as pickle

  # If the input is a string then open and read from that file
  if isinstance(infile, str):
    with open(infile, 'rb') as infile:
      return pickle.load(infile)

  # Otherwise assume the input is a file and read from it
  else:
    return pickle.load(infile)
