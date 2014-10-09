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


class Script(object):
  '''Class to run script.'''

  def __init__(self, imagesets, reflections, crystals=None):
    '''Setup the script.'''

    # Filename data
    self.imagesets = imagesets
    self.reflections = reflections
    self.crystals = crystals

  def __call__(self):
    '''Run the script.'''
    from dials.array_family import flex # import dependency

    self.view()

  def view(self):
    from dials.util.spotfinder_wrap import spot_wrapper
    spot_wrapper(working_phil=None).display(
      imagesets=self.imagesets, reflections=self.reflections,
      crystals=self.crystals)

if __name__ == '__main__':
  import sys

  from dials.util.options import OptionParser
  from dials.util.options import flatten_datablocks
  from dials.util.options import flatten_experiments
  from dials.util.options import flatten_reflections
  import libtbx.load_env
  usage_message = """
    %s datablock.json reflections.pickle
  """ %libtbx.env.dispatcher_name
  parser = OptionParser(
    usage=usage_message,
    read_datablocks=True,
    read_experiments=True,
    read_reflections=True,
    read_datablocks_from_images=True)
  params, options = parser.parse_args(show_diff_phil=True)
  datablocks = flatten_datablocks(params.input.datablock)
  experiments = flatten_experiments(params.input.experiments)
  reflections = flatten_reflections(params.input.reflections)

  if len(datablocks) == 0 and len(experiments) == 0:
    parser.print_help()
    exit(0)

  if len(datablocks) > 0:
    assert(len(datablocks) == 1)
    imagesets = datablocks[0].extract_imagesets()
    crystals = None
  elif len(experiments.imagesets()) > 0:
    assert(len(experiments.imagesets()) == 1)
    imagesets = experiments.imagesets()
    crystals = experiments.crystals()
  else:
    raise RuntimeError("No imageset could be constructed")

  runner = Script(
      reflections=reflections,
      imagesets=imagesets,
      crystals=crystals)

  # Run the script
  runner()
