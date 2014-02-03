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

from __future__ import division

# Import to give access from here
from dxtbx.serialize.dump import imageset as sweep
from dxtbx.serialize.dump import imageset_to_string as sweep_to_string
from dxtbx.serialize.dump import datablock

# FIXME These are only temporary - they should be removed ASAP
from cctbx.crystal.crystal_model.serialize import crystal_to_string
from cctbx.crystal.crystal_model.serialize import dump_crystal as crystal


def reflections(obj, outfile):
  ''' Dump the given object to file

  Params:
      obj The reflection list to dump
      outfile The output file name or file object

  '''
  import cPickle as pickle

  if isinstance(outfile, str):
    with open(outfile, 'wb') as outfile:
      pickle.dump(obj, outfile, pickle.HIGHEST_PROTOCOL)

  # Otherwise assume the input is a file and write to it
  else:
    pickle.dump(obj, outfile, pickle.HIGHEST_PROTOCOL)

def reference(obj, outfile):
  ''' Dump the given object to file

  Params:
      obj The reference list to dump
      outfile The output file name or file object

  '''
  import cPickle as pickle

  if isinstance(outfile, str):
    with open(outfile, 'wb') as outfile:
      pickle.dump(obj, outfile, pickle.HIGHEST_PROTOCOL)

  # Otherwise assume the input is a file and write to it
  else:
    pickle.dump(obj, outfile, pickle.HIGHEST_PROTOCOL)

def experiment_list(obj, outfile):
  ''' Dump an experiment list. '''
  from dials.model.experiment.experiment_list import ExperimentListDumper
  dumper = ExperimentListDumper(obj)
  dumper.as_file(outfile)
