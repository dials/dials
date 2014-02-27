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
from iotbx.phil import parse

phil_scope = parse('''

spotfinder
  .help = "Parameters used in the spot finding algorithm."
{
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

  threshold
    .help = "Parameters used in the spot finding threshold strategy."
  {
    algorithm = *xds unimodal lui
      .help = "The threshold strategy."
      .type = choice

    kernel_size = 3 3
      .help = "The size of the local area around the spot in which"
              "to calculate the mean and variance. The kernel is"
              "given as a box of size (2 * nx + 1, 2 * ny + 1) centred"
              "at the pixel."
              ""
              "Used by: xds."
      .type = ints(size=2)

    sigma_background = 6
      .help = "The number of standard deviations of the coefficient of"
              "variation (variance / mean) in the local area below"
              "which the pixel will be classified as background."
              ""
              "Used by: xds."
      .type = float

    sigma_strong = 3
      .help = "The number of standard deviations above the mean in the"
              "local area above which the pixel will be classified as"
              "strong."
              ""
              "Used by: xds."
      .type = float

    min_local=2
      .help = "The number of pixels in the local area of each pixel needed"
              "to do the thresholding. Setting to 0 or less means that all"
              "the pixels under the kernel are needed. The minimum allowable"
              "number is 2"
      .type = int

    times = 5
      .help = "How many times to run the smoothing algorithm."
              ""
              "Used by: lui"
      .type = int

    block_size = 5 12
      .help = "How many (nx x ny) blocks per image."
              ""
              "Used by: lui"
      .type = ints(size=2)

    shift = 10
      .help = "The shift value (in counts) added to the threshold."
              ""
              "Used by: lui"
      .type = int

    dimensions = *2d 3d
      .help = "Use the 2D or 3D smoothing algorithm."
              ""
              "Used by: lui"
      .type = choice
  }

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

''')
