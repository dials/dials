#!/usr/bin/env python
#
# spotfinding.py
#
#  Copyright (C) 2013 Diamond Light Source
#
#  Author: James Parkhurst
#
#  This code is distributed under the BSD license, a copy of which is
#  included in the root directory of this package.
from __future__ import division

def generate_phil_scope():
  from iotbx.phil import parse
  from dials.interfaces import SpotFinderThresholdIface

  phil_scope = parse('''

  spotfinder
    .help = "Parameters used in the spot finding algorithm."
  {
    include scope dials.data.lookup.phil_scope

    scan_range = None
      .help = "The range of images to use in finding spots. Number of arguments"
              "must be a factor of two. Specifying \"0 0\" will use all images"
              "by default. The given range follows C conventions"
              "(e.g. j0 <= j < j1)."
              "For sweeps the scan range is interpreted as the literal scan"
              "range. Whereas for imagesets the scan range is interpreted as"
              "the image number in the imageset"
      .type = ints(size=2)
      .multiple = True

    filter
      .help = "Parameters used in the spot finding filter strategy."

    {
      min_spot_size = 6
        .help = "The minimum number of contiguous pixels for a spot"
                "to be accepted by the filtering algorithm."
                ""
                "Used by: xds."
        .type = int(value_min=0)

      max_separation = 2
        .help = "The maximum peak-to-centroid separation (in pixels)"
                "for a spot to be accepted by the filtering algorithm."
                ""
                "Used by: xds."
        .type = float(value_min=0)

      d_min = None
        .help = "The high resolution limit in Angstrom for a spot to be"
                "accepted by the filtering algorithm."
        .type = float(value_min=0)

      d_max = None
        .help = "The low resolution limit in Angstrom for a spot to be"
                "accepted by the filtering algorithm."
        .type = float(value_min=0)

      ice_rings {
        filter = False
          .type = bool
        unit_cell = 4.498,4.498,7.338,90,90,120
          .type = unit_cell
          .help = "The unit cell to generate d_spacings for powder rings."
        space_group = 194
          .type = space_group
          .help = "The space group used to generate d_spacings for powder rings."
      }

      untrusted_polygon = None
        .multiple = True
        .type = ints(value_min=0)

      #untrusted_ellipse = None
      #  .multiple = True
      #  .type = ints(size=4, value_min=0)

      #untrusted_rectangle = None
      #  .multiple = True
      #  .type = ints(size=4, value_min=0)
    }
    save_shoeboxes = True
      .type = bool
      .help = "Save the raw pixel values inside the reflection shoeboxes."
    image_viewer = False
      .help = "Display the results of spot finding in an interactive image"
              "viewer."
      .type = bool
  }

  ''', process_includes=True)

  main_scope = phil_scope.get_without_substitution("spotfinder")
  assert(len(main_scope) == 1)
  main_scope = main_scope[0]
  main_scope.adopt_scope(SpotFinderThresholdIface.phil_scope())
  return phil_scope

phil_scope = generate_phil_scope()
