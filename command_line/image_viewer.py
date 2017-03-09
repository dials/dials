#!/usr/bin/env python
#
# image_viewer.py
#
#  Copyright (C) 2013 Diamond Light Source
#
#  Author: Richard Gildea
#
#  This code is distributed under the BSD license, a copy of which is
#  included in the root directory of this package.

from __future__ import absolute_import, division
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
show_indexed = False
  .type = bool
show_integrated = False
  .type = bool
show_mask = False
  .type = bool
show_mask2 = False
  .type = bool
show_basis_vectors = True
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
d_min = None
  .type = float(value_min=0)
mask = None
  .type = str
  .help = path to mask pickle file

masking {
  include scope dials.util.masking.phil_scope
}

output {
  mask = mask.pickle
    .type = str
    .help = "Name of output mask file"
}

predict_reflections = False
  .type = bool
  .help = Predict reflections if no reflections provided in input

include scope dials.algorithms.profile_model.factory.phil_scope
include scope dials.algorithms.spot_prediction.reflection_predictor.phil_scope

""", process_includes=True)

class Script(object):
  '''Class to run script.'''

  def __init__(self, params, datablock, experiments, reflections):
    '''Setup the script.'''

    # Filename data
    self.params = params
    self.datablock = datablock
    self.experiments = experiments
    self.reflections = reflections
    self.wrapper = None

  def __call__(self):
    '''Run the script.'''
    from dials.array_family import flex # import dependency

    self.view()

  def view(self):
    from dials.util.spotfinder_wrap import spot_wrapper
    self.wrapper = spot_wrapper(params=self.params)
    self.wrapper.display(
      datablock=self.datablock,
      experiments=self.experiments,
      reflections=self.reflections)

if __name__ == '__main__':
  import wx # It is unclear why, but it is crucial that wx
            # is imported before the parser is run.
            # Otherwise viewer will crash when run with
            # .cbf image as parameter on linux with wxPython>=3
            # The problem can be traced to
            # dxtbx/format/FormatCBFFull.py:49
            #  ''' from iotbx.detectors.cbf import CBFImage '''
            # and the wx import must happen before that import.

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
    assert len(datablocks) == 1
    datablock = datablocks[0]
  else:
    datablock = None

  if params.mask is not None:
    from libtbx import easy_pickle
    params.mask = easy_pickle.load(params.mask)

  runner = Script(
    params=params,
    reflections=reflections,
    datablock=datablock,
    experiments=experiments
  )

  # Run the script
  runner()
