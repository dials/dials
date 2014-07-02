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

  def __init__(self, imagesets, reflections, crystals=None):
    '''Setup the script.'''

    # Filename data
    self.imagesets = imagesets
    self.reflections = reflections
    self.crystals = crystals

  def __call__(self):
    '''Run the script.'''
    from dials.model.data import ReflectionList # import dependency

    self.view()

  def view(self):
    from dials.util.spotfinder_wrap import spot_wrapper
    spot_wrapper(working_phil=None).display(
      imagesets=self.imagesets, reflections=self.reflections,
      crystals=self.crystals)

if __name__ == '__main__':
  import sys
  from dials.util.command_line import Importer
  args = sys.argv[1:]
  if len(args) == 0:
    from libtbx.utils import Usage
    import libtbx.load_env
    usage_message = """\
%s datablock.json reflections.pickle""" %libtbx.env.dispatcher_name
    raise Usage(usage_message)
  importer = Importer(args, check_format=True)
  crystals = None
  if importer.datablocks is not None and len(importer.datablocks) == 1:
    imagesets = importer.datablocks[0].extract_imagesets()
  elif importer.datablocks is not None and len(importer.datablocks) > 1:
    raise RuntimeError("Only one DataBlock can be processed at a time")
  elif len(importer.experiments.imagesets()) > 0:
    imagesets = importer.experiments.imagesets()[:1]
    crystals = importer.experiments.crystals()
  else:
    raise RuntimeError("No imageset could be constructed")
  assert len(importer.unhandled_arguments) == 0

  runner = ScriptRunner(
      reflections=importer.reflections,
      imagesets=imagesets,
      crystals=crystals)

  # Run the script
  runner()
