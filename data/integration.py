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

def generate_phil_scope():

  from libtbx.phil import parse
  from dials.interfaces import CentroidIface
  from dials.interfaces import BackgroundIface
  from dials.interfaces import IntensityIface

  phil_scope = parse('''

  integration
    .help = "Configure the integration algorithm."
  {
    include scope dials.data.lookup.phil_scope

    shoebox
      .help = "Parameters used in the integration shoebox"
    {
      block_size = 10
        .help = "The block size in rotation angle (degrees)."
        .type = float

      n_sigma = 3
        .help = "The number of standard deviations of the beam divergence and the"
                "mosaicity to use for the bounding box size."
        .type = float

      sigma_b = None
        .help = "The E.S.D. of the beam divergence"
        .type = float

      sigma_m = None
        .help = "The E.S.D. of the reflecting range"
        .type = float
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

  ''', process_includes=True)

  main_scope = phil_scope.get_without_substitution("integration")
  assert(len(main_scope) == 1)
  main_scope = main_scope[0]
  main_scope.adopt_scope(CentroidIface.phil_scope())
  main_scope.adopt_scope(BackgroundIface.phil_scope())
  main_scope.adopt_scope(IntensityIface.phil_scope())

  return phil_scope

phil_scope = generate_phil_scope()
