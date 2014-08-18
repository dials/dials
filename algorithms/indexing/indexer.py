#!/usr/bin/env python
# -*- mode: python; coding: utf-8; indent-tabs-mode: nil; python-indent: 2 -*-
#
# dials.algorithms.indexing.indexer.py
#
#  Copyright (C) 2013 Diamond Light Source
#
#  Author: Graeme Winter, Nick Sauter, Richard Gildea
#
#  This code is distributed under the BSD license, a copy of which is
#  included in the root directory of this package.

from __future__ import division
import math
from scitbx import matrix

from rstbx.array_family import flex
import copy


# we need these things

def discover_better_experimental_model(spot_positions, detector, beam,
                                       goniometer, scan, params):
  '''Given an attempt at indexing derive a more likely model for the
  experimental geometry.'''

  # Spot_positions: Centroid positions for spotfinder spots, in pixels
  # Return value: Corrected for parallax, converted to mm

  from dials.algorithms.indexing.indexer2 import indexer_base
  spots_mm = indexer_base.map_spots_pixel_to_mm_rad(
    spots=spot_positions, detector=detector, scan=scan)

  # derive a max_cell from mm spots
  # derive a grid sampling from spots

  from rstbx.indexing_api.lattice import DPS_primitive_lattice
  # max_cell: max possible cell in Angstroms; set to None, determine from data
  # recommended_grid_sampling_rad: grid sampling in radians; guess for now

  DPS = DPS_primitive_lattice(max_cell = None,
                              recommended_grid_sampling_rad = None,
                              horizon_phil = params)
  from scitbx import matrix
  DPS.S0_vector = matrix.col(beam.get_s0())
  DPS.inv_wave = 1./beam.get_wavelength()
  if goniometer is None:
    DPS.axis = matrix.col((1,0,0))
  else:
    DPS.axis = matrix.col(goniometer.get_rotation_axis())
  DPS.set_detector(detector)

  # transform input into what Nick needs
  # i.e., construct a flex.vec3 double consisting of mm spots, phi in degrees

  data = flex.vec3_double()
  for spot in spots_mm:
    data.append((spot['xyzobs.mm.value'][0],
                 spot['xyzobs.mm.value'][1],
                 spot['xyzobs.mm.value'][2]*180./math.pi))

  #from matplotlib import pyplot as plt
  #plt.plot([spot.centroid_position[0] for spot in spots_mm] , [spot.centroid_position[1] for spot in spots_mm], 'ro')
  #plt.show()

  DPS.index(raw_spot_input = data, panel_addresses = flex.int([s['panel'] for s in spots_mm]))

  # for development, we want an exhaustive plot of beam probability map:
  params.indexing.plot_search_scope = False

  # perform calculation
  if params.indexing.improve_local_scope=="origin_offset":
    new_detector = DPS.optimize_origin_offset_local_scope()
    return new_detector, beam
  elif params.indexing.improve_local_scope=="S0_vector":
    new_S0_vector = DPS.optimize_S0_local_scope()
    import copy
    new_beam = copy.copy(beam)
    new_beam.set_s0(new_S0_vector)
    return detector, new_beam

  #NKS TO DO: implement rigorous scope instead of local scope.
