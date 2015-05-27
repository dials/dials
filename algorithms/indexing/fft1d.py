#!/usr/bin/env python
# -*- mode: python; coding: utf-8; indent-tabs-mode: nil; python-indent: 2 -*-
#
# dials.algorithms.indexing.fft1d.py
#
#  Copyright (C) 2014 Diamond Light Source
#
#  Author: Richard Gildea
#
#  This code is distributed under the BSD license, a copy of which is
#  included in the root directory of this package.

from __future__ import division

import math
from scitbx.array_family import flex
from dials.algorithms.indexing.indexer import indexer_base
from dxtbx.model.experiment.experiment_list import Experiment, ExperimentList

class indexer_fft1d(indexer_base):

  def __init__(self, reflections, sweep, params):
    super(indexer_fft1d, self).__init__(reflections, sweep, params)

  def find_lattices(self):
    self.d_min = self.params.refinement_protocol.d_min_start

    from rstbx.phil.phil_preferences import indexing_api_defs
    import iotbx.phil
    hardcoded_phil = iotbx.phil.parse(
      input_string=indexing_api_defs).extract()

    reflections = self.reflections.select(
      (self.reflections['id'] == -1) &
      (1/self.reflections['rlp'].norms() > self.d_min))
    solutions = candidate_basis_vectors_fft1d(
      reflections, self.detector, self.beam,
      self.goniometer, self.imagesets[0].get_scan(), hardcoded_phil,
      max_cell=self.params.max_cell)
    self.candidate_basis_vectors = solutions[0]
    if self.params.debug:
      self.debug_show_candidate_basis_vectors()
    if self.params.debug_plots:
      self.debug_plot_candidate_basis_vectors()
    self.candidate_crystal_models = self.find_candidate_orientation_matrices(
      self.candidate_basis_vectors,
      max_combinations=self.params.basis_vector_combinations.max_try)
    crystal_model, n_indexed = self.choose_best_orientation_matrix(
      self.candidate_crystal_models)
    if crystal_model is not None:
      crystal_models = [crystal_model]
    else:
      crystal_models = []
    experiments = ExperimentList()
    for cm in crystal_models:
      experiments.append(Experiment(beam=self.beam,
                                    detector=self.detector,
                                    goniometer=self.goniometer,
                                    scan=self.imagesets[0].get_scan(),
                                    crystal=cm))
    return experiments

def candidate_basis_vectors_fft1d(spot_positions, detector, beam,
                                  goniometer, scan, params,
                                  max_cell=None):

  # Spot_positions: Centroid positions for spotfinder spots, in pixels
  # Return value: Corrected for parallax, converted to mm

  spots_mm = indexer_base.map_spots_pixel_to_mm_rad(
    spots=spot_positions, detector=detector, scan=scan)

  if len(detector) > 1:
    panel_ids = [spot['panel'] for spot in spot_positions]
  else:
    panel_ids = None


  # derive a max_cell from mm spots
  # derive a grid sampling from spots

  from rstbx.indexing_api.lattice import DPS_primitive_lattice
  # max_cell: max possible cell in Angstroms; set to None, determine from data
  # recommended_grid_sampling_rad: grid sampling in radians; guess for now

  DPS = DPS_primitive_lattice(max_cell=max_cell,
                              recommended_grid_sampling_rad = None,
                              horizon_phil = params)
  from scitbx import matrix
  DPS.S0_vector = matrix.col(beam.get_s0())
  DPS.inv_wave = 1./beam.get_wavelength()

  # transform input into what Nick needs
  # i.e., construct a flex.vec3 double consisting of mm spots, phi in degrees

  DPS.index(reciprocal_space_vectors=spots_mm['rlp'])
  solutions = DPS.getSolutions()
  return [matrix.col(s.bvec()) for s in solutions],DPS.getXyzData()
