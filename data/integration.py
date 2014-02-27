#!/usr/bin/env python
#
# integration.py
#
#  Copyright (C) 2013 Diamond Light Source
#
#  Author: James Parkhurst
#
#  This code is distributed under the BSD license, a copy of which is
#  included in the root directory of this package.
from __future__ import division
from libtbx.phil import parse

phil_scope = parse('''

integration
  .help = "Configure the integration algorithm."
{
  algorithm = sum2d *sum3d mosflm fitrs
    .help = "The integration algorithm"
    .type = choice

  fitrs
    .help = "Parameters for reciprocal space profile fitting."
  {
    profile
      .help = "Parameters for profile fitting"
    {
      grid_size = 5
        .help = "The size of the reciprocal space grid for each reflection."
                "The size is the same in each dimensions"
        .type = int

      reference_frame_interval = 10
        .help = "The oscillation at which to learn new reference profiles"
        .type = int

      reference_signal_threshold = 0.02
        .help = "The threshold to use in reference profile"
        .type = float
    }
  }

  mosflm
    .help = "Parameters for mosflm profile fitting"
  {
    nblocks = 4
      .help = "number of block per coordinate"
      .type = int
  }

  filter
    .help = "Parameters for filtering reflections"
  {
    by_bbox = False
      .help = "Filter the reflections by the volume of the bounding box."
              "A threshold value is chosen from a histogram of the volumes."
              "Reflections with bounding box volume above the threshold value"
              "are not used in intergration."
      .type = bool

    by_zeta = 0.05
      .help = "Filter the reflections by the value of zeta. A value of less"
              "than or equal to zero indicates that this will not be used. A"
              "positive value is used as the minimum permissable value."
      .type = float

    by_xds_small_angle = False
      .help = "Filter the reflections by whether the XDS small angle"
              "approximation holds for the reflection."
      .type = bool

    by_xds_angle = False
      .help = "Filter the reflections by whether the geometry of the XDS"
              "transform allows a reflection to be transformed."
      .type = bool
  }
}

''')
