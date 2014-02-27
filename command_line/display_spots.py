#!/usr/bin/env python
#
# display_spots.py
#
#  Copyright (C) 2013 Diamond Light Source
#
#  Author: Richard Gildea
#
#  This code is distributed under the BSD license, a copy of which is
#  included in the root directory of this package.

from __future__ import division
# LIBTBX_PRE_DISPATCHER_INCLUDE_SH export PHENIX_GUI_ENVIRONMENT=1
# LIBTBX_PRE_DISPATCHER_INCLUDE_SH export BOOST_ADAPTBX_FPE_DEFAULT=1


class ScriptRunner(object):
  '''Class to run script.'''

  def __init__(self, sweep_filenames, reflections):
    '''Setup the script.'''

    # Filename data
    self.sweep_filenames = sweep_filenames
    self.reflections = reflections

  def __call__(self):
    '''Run the script.'''
    import cPickle as pickle
    from dials.model.data import ReflectionList # import dependency
    from dials.util.command_line import Command

    self.view()

  def view(self):
    from dials.util.spotfinder_wrap import spot_wrapper
    spot_wrapper(working_phil=None).display(
        sweep_filenames=self.sweep_filenames, reflections=self.reflections)

if __name__ == '__main__':
  import sys
  from optparse import OptionParser

  # Specify the command line options
  usage  = "usage: %prog [options] " \
           "/path/to/reflections.pickle " \
           "/path/to/image/files "

  # Create an option parser
  parser = OptionParser(usage)

  # Parse the arguments
  options, args = parser.parse_args()

  # Print help if no arguments specified, otherwise call function
  if len(args) < 1:
    parser.print_help()
  else:
    # Initialise the script runner
    from dials.util.command_line import Importer

    args = sys.argv[1:]
  importer = Importer(args, check_format=False)
  if importer.datablocks is not None and len(importer.datablocks) == 1:
    imagesets = importer.datablocks[0].extract_imagesets()
  elif importer.datablocks is not None and len(importer.datablocks) > 1:
    raise RuntimeError("Only one DataBlock can be processed at a time")
  elif len(importer.experiments.imagesets()) == 1:
    imagesets = importer.experiments.imagesets()
  else:
    raise RuntimeError("No imageset could be constructed")
  paths = []
  for imageset in imagesets:
    paths.extend(imageset.paths())
  assert len(importer.unhandled_arguments) == 0

  runner = ScriptRunner(
      reflections=importer.reflections,
      sweep_filenames=paths)

  # Run the script
  runner()
