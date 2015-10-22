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


import iotbx.phil

help_message = '''

This program can be used for viewing diffraction images, optionally overlayed
with the results of spot finding, indexing or integration.

Examples::

  dials.image_viewer image.cbf

  dials.image_viewer datablock.json

  dials.image_viewer datablock.json strong.pickle

  dials.image_viewer datablock.json integrated.pickle

  dials.image_viewer experiments.json

'''

phil_scope = iotbx.phil.parse("""\
image_viewer {
  brightness = 100
    .type = int
  color_scheme = *grayscale rainbow heatmap invert
    .type = choice
  show_beam_center = True
    .type = bool
  show_resolution_rings = False
    .type = bool
  show_ice_rings = False
    .type = bool
  show_ctr_mass = True
    .type = bool
  show_max_pix = True
    .type = bool
  show_all_pix = True
    .type = bool
  show_shoebox = True
    .type = bool
  show_predictions = True
    .type = bool
  show_miller_indices = False
    .type = bool
  display = *image mean variance dispersion sigma_b \
            sigma_s threshold global_threshold
    .type = choice
  nsigma_b = 6
    .type = float(value_min=0)
  nsigma_s = 3
    .type = float(value_min=0)
  global_threshold = 0
    .type = float(value_min=0)
  kernel_size = 3,3
    .type = ints(size=2, value_min=1)
  min_local = 2
    .type = int
  gain = 1
    .type = float(value_min=0)
  sum_images = 1
    .type = int(value_min=1)
    .expert_level = 2
  untrusted_polygon = None
    .multiple = True
    .type = ints(value_min=0)
}
""")

class Script(object):
  '''Class to run script.'''

  def __init__(self, params, imagesets, reflections, crystals=None):
    '''Setup the script.'''

    # Filename data
    self.params = params
    self.imagesets = imagesets
    self.reflections = reflections
    self.crystals = crystals

  def __call__(self):
    '''Run the script.'''
    from dials.array_family import flex # import dependency

    self.view()

  def view(self):
    from dials.util.spotfinder_wrap import spot_wrapper
    spot_wrapper(params=self.params).display(
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
    %s datablock.json [reflections.pickle]
  """ %libtbx.env.dispatcher_name
  parser = OptionParser(
    usage=usage_message,
    phil=phil_scope,
    read_datablocks=True,
    read_experiments=True,
    read_reflections=True,
    read_datablocks_from_images=True,
    epilog=help_message)
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
    params=params.image_viewer,
    reflections=reflections,
    imagesets=imagesets,
    crystals=crystals)

  # Run the script
  runner()
