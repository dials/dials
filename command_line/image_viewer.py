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

from __future__ import absolute_import, division, print_function
import sys

import dials.util.banner  # noqa: F401 - Importing means that it prints
import iotbx.phil

# LIBTBX_PRE_DISPATCHER_INCLUDE_SH export PHENIX_GUI_ENVIRONMENT=1
# LIBTBX_PRE_DISPATCHER_INCLUDE_SH export BOOST_ADAPTBX_FPE_DEFAULT=1


help_message = """

This program can be used for viewing diffraction images, optionally overlayed
with the results of spot finding, indexing or integration.

Examples::

  dials.image_viewer image.cbf

  dials.image_viewer models.expt

  dials.image_viewer models.expt strong.refl

  dials.image_viewer models.expt integrated.refl

"""

phil_scope = iotbx.phil.parse(
    """\
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
  .help = "Set gain for the thresholding algorithm. This does not override the"
          "detector's panel gain, but acts as a multiplier for it."
sum_images = 1
  .type = int(value_min=1)
  .expert_level = 2
d_min = None
  .type = float(value_min=0)
mask = None
  .type = path
  .help = path to mask pickle file

include scope rstbx.phil.phil_preferences.iotbx_defs_viewer_detail

masking {
  include scope dials.util.masking.phil_scope
}

output {
  mask = pixels.mask
    .type = path
    .help = "Name of output mask file"

  mask_params = mask.phil
    .type = path
    .help = "Name of output mask parameter file"
}

predict_reflections = False
  .type = bool
  .help = Predict reflections if no reflections provided in input

include scope dials.algorithms.profile_model.factory.phil_scope
include scope dials.algorithms.spot_prediction.reflection_predictor.phil_scope

""",
    process_includes=True,
)


class Script(object):
    """Class to run script."""

    def __init__(self, params, experiments, reflections):
        """Setup the script."""

        # Filename data
        self.params = params
        self.experiments = experiments
        self.reflections = reflections
        self.wrapper = None

    def __call__(self):
        """Run the script."""
        from dials.array_family import flex  # noqa

        self.view()

    def view(self):
        from dials.util.image_viewer.spotfinder_wrap import spot_wrapper

        self.wrapper = spot_wrapper(params=self.params)
        self.wrapper.display(experiments=self.experiments, reflections=self.reflections)


if __name__ == "__main__":
    import wx  # It is unclear why, but it is crucial that wx

    # is imported before the parser is run.
    # Otherwise viewer will crash when run with
    # .cbf image as parameter on linux with wxPython>=3
    # The problem can be traced to
    # dxtbx/format/FormatCBFFull.py:49
    #  ''' from iotbx.detectors.cbf import CBFImage '''
    # and the wx import must happen before that import.
    WX3 = wx.VERSION[0] == 3
    if not WX3:
        # HACK: Monkeypatch this renamed function so we can trick wxtbx's IntCtrl
        #       without having to alter the package
        wx.SystemSettings_GetColour = wx.SystemSettings.GetColour

    from dials.util.options import OptionParser
    from dials.util.options import flatten_reflections
    from dials.util.options import flatten_experiments

    import libtbx.load_env

    usage_message = (
        """
    %s models.expt [observations.refl]
  """
        % libtbx.env.dispatcher_name
    )
    parser = OptionParser(
        usage=usage_message,
        phil=phil_scope,
        read_experiments=True,
        read_reflections=True,
        read_experiments_from_images=True,
        epilog=help_message,
    )
    params, options = parser.parse_args(show_diff_phil=True)
    experiments = [x.data for x in params.input.experiments]
    reflections = flatten_reflections(params.input.reflections)

    if len(experiments) == 0:
        parser.print_help()
        exit(0)

    flat_expts = flatten_experiments(params.input.experiments)
    if not all(e.detector for e in flat_expts):
        sys.exit("Error: experiment has no detector")
    if not all(e.beam for e in flat_expts):
        sys.exit("Error: experiment has no beam")

    if params.mask is not None:
        from libtbx import easy_pickle

        params.mask = easy_pickle.load(params.mask)

    runner = Script(params=params, reflections=reflections, experiments=experiments)

    # Run the script
    runner()
