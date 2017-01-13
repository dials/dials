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

  def __init__(self, reflections, imagesets, params):
    super(indexer_fft1d, self).__init__(reflections, imagesets, params)

  def find_candidate_basis_vectors(self):
    self.d_min = self.params.refinement_protocol.d_min_start

    from rstbx.phil.phil_preferences import indexing_api_defs
    import iotbx.phil
    hardcoded_phil = iotbx.phil.parse(
      input_string=indexing_api_defs).extract()

    sel = (self.reflections['id'] == -1)
    if self.d_min is not None:
      sel &= (1/self.reflections['rlp'].norms() > self.d_min)
    reflections = self.reflections.select(sel)
    solutions = candidate_basis_vectors_fft1d(
      reflections['rlp'], hardcoded_phil, max_cell=self.params.max_cell)
    self.candidate_basis_vectors = solutions[0]
    self.debug_show_candidate_basis_vectors()
    if self.params.debug_plots:
      self.debug_plot_candidate_basis_vectors()

    return self.candidate_basis_vectors

  def find_lattices(self):
    self.find_candidate_basis_vectors()
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
      for imageset in self.imagesets:
        experiments.append(Experiment(imageset=imageset,
                                      beam=imageset.get_beam(),
                                      detector=imageset.get_detector(),
                                      goniometer=imageset.get_goniometer(),
                                      scan=imageset.get_scan(),
                                      crystal=cm))
    return experiments

def candidate_basis_vectors_fft1d(reciprocal_lattice_vectors, params,
                                  max_cell=None):

  # Spot_positions: Centroid positions for spotfinder spots, in pixels
  # Return value: Corrected for parallax, converted to mm

  # derive a max_cell from mm spots
  # derive a grid sampling from spots

  from rstbx.indexing_api.lattice import DPS_primitive_lattice
  # max_cell: max possible cell in Angstroms; set to None, determine from data
  # recommended_grid_sampling_rad: grid sampling in radians; guess for now

  DPS = DPS_primitive_lattice(max_cell=max_cell,
                              recommended_grid_sampling_rad = None,
                              horizon_phil = params)
  from scitbx import matrix
  #DPS.S0_vector = matrix.col(beam.get_s0())
  #DPS.inv_wave = 1./beam.get_wavelength()

  # transform input into what Nick needs
  # i.e., construct a flex.vec3 double consisting of mm spots, phi in degrees

  DPS.index(reciprocal_space_vectors=reciprocal_lattice_vectors)
  solutions = DPS.getSolutions()
  return [matrix.col(s.bvec()) for s in solutions],DPS.getXyzData()
