#!/usr/bin/env python
#
# dials.model.serialize.dump.py
#
#  Copyright (C) 2013 Diamond Light Source
#
#  Author: James Parkhurst
#
#  This code is distributed under the BSD license, a copy of which is
#  included in the root directory of this package.

from __future__ import absolute_import, division

# Import to give access from here
from dxtbx.serialize.dump import imageset as sweep # implicit import
from dxtbx.serialize.dump import imageset_to_string as sweep_to_string # implicit import
from dxtbx.serialize.dump import datablock # implicit import


def reflections(obj, outfile):
  '''
  Dump the given object to file

  :param obj: The reflection list to dump
  :param outfile: The output file name or file object

  '''
  import cPickle as pickle

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
  import cPickle as pickle

  if isinstance(outfile, str):
    with open(outfile, 'wb') as outfile:
      pickle.dump(obj, outfile, pickle.HIGHEST_PROTOCOL)

  # Otherwise assume the input is a file and write to it
  else:
    pickle.dump(obj, outfile, pickle.HIGHEST_PROTOCOL)
