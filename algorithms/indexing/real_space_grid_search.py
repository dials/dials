#!/usr/bin/env python
# -*- mode: python; coding: utf-8; indent-tabs-mode: nil; python-indent: 2 -*-
#
# dials.algorithms.indexing.real_space_grid_search.py
#
#  Copyright (C) 2014 Diamond Light Source
#
#  Author: Richard Gildea
#
#  This code is distributed under the BSD license, a copy of which is
#  included in the root directory of this package.

from __future__ import absolute_import, division
import math
import logging
logger = logging.getLogger(__name__)

from scitbx import matrix
from scitbx.array_family import flex
from dials.algorithms.indexing.indexer import \
     indexer_base, optimise_basis_vectors
from dials.algorithms.indexing.indexer import \
     is_approximate_integer_multiple
from dxtbx.model.experiment_list import Experiment, ExperimentList


class indexer_real_space_grid_search(indexer_base):

  def __init__(self, reflections, imagesets, params):
    super(indexer_real_space_grid_search, self).__init__(
      reflections, imagesets, params)

  def find_lattices(self):
    self.real_space_grid_search()
    crystal_models = self.candidate_crystal_models
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

  def real_space_grid_search(self):
    d_min = self.params.refinement_protocol.d_min_start

    sel = (self.reflections['id'] == -1)
    if d_min is not None:
      sel &= (1/self.reflections['rlp'].norms() > d_min)
    reciprocal_lattice_points = self.reflections['rlp'].select(sel)

    logger.info("Indexing from %i reflections" %len(reciprocal_lattice_points))

    def compute_functional(vector):
      two_pi_S_dot_v = 2 * math.pi * reciprocal_lattice_points.dot(vector)
      return flex.sum(flex.cos(two_pi_S_dot_v))

    from rstbx.array_family import flex
    from rstbx.dps_core import SimpleSamplerTool
    assert self.target_symmetry_primitive is not None
    assert self.target_symmetry_primitive.unit_cell() is not None
    SST = SimpleSamplerTool(
      self.params.real_space_grid_search.characteristic_grid)
    SST.construct_hemisphere_grid(SST.incr)
    cell_dimensions = self.target_symmetry_primitive.unit_cell().parameters()[:3]
    unique_cell_dimensions = set(cell_dimensions)
    logger.info(
      "Number of search vectors: %i" %(len(SST.angles) * len(unique_cell_dimensions)))
    vectors = flex.vec3_double()
    function_values = flex.double()
    for i, direction in enumerate(SST.angles):
      for l in unique_cell_dimensions:
        v = matrix.col(direction.dvec) * l
        f = compute_functional(v.elems)
        vectors.append(v.elems)
        function_values.append(f)

    perm = flex.sort_permutation(function_values, reverse=True)
    vectors = vectors.select(perm)
    function_values = function_values.select(perm)

    unique_vectors = []
    i = 0
    while len(unique_vectors) < 30:
      v = matrix.col(vectors[i])
      is_unique = True
      if i > 0:
        for v_u in unique_vectors:
          if v.length() < v_u.length():
            if is_approximate_integer_multiple(v, v_u):
              is_unique = False
              break
          elif is_approximate_integer_multiple(v_u, v):
            is_unique = False
            break
      if is_unique:
        unique_vectors.append(v)
      i += 1

    for i in range(30):
      v = matrix.col(vectors[i])
      logger.debug("%s %s %s" %(str(v.elems), str(v.length()), str(function_values[i])))

    basis_vectors = [v.elems for v in unique_vectors]
    self.candidate_basis_vectors = basis_vectors

    if self.params.optimise_initial_basis_vectors:
      optimised_basis_vectors = optimise_basis_vectors(
        reciprocal_lattice_points, basis_vectors)
      optimised_function_values = flex.double([
        compute_functional(v) for v in optimised_basis_vectors])

      perm = flex.sort_permutation(optimised_function_values, reverse=True)
      optimised_basis_vectors = optimised_basis_vectors.select(perm)
      optimised_function_values = optimised_function_values.select(perm)

      unique_vectors = [matrix.col(v) for v in optimised_basis_vectors]

    logger.info("Number of unique vectors: %i" %len(unique_vectors))

    for i in range(len(unique_vectors)):
      logger.debug("%s %s %s" %(
        str(compute_functional(unique_vectors[i].elems)),
        str(unique_vectors[i].length()),
        str(unique_vectors[i].elems)))

    crystal_models = []
    self.candidate_basis_vectors = unique_vectors
    self.debug_show_candidate_basis_vectors()
    if self.params.debug_plots:
      self.debug_plot_candidate_basis_vectors()
    candidate_orientation_matrices \
      = self.find_candidate_orientation_matrices(unique_vectors)
    crystal_model, n_indexed = self.choose_best_orientation_matrix(
      candidate_orientation_matrices)
    if crystal_model is not None:
      crystal_models = [crystal_model]
    else:
      crystal_models = []

    #assert len(crystal_models) > 0

    candidate_orientation_matrices = crystal_models

    #for i in range(len(candidate_orientation_matrices)):
      #if self.target_symmetry_primitive is not None:
        ##print "symmetrizing model"
        ##self.target_symmetry_primitive.show_summary()
        #symmetrized_model = self.apply_symmetry(
          #candidate_orientation_matrices[i], self.target_symmetry_primitive)
        #candidate_orientation_matrices[i] = symmetrized_model

    self.candidate_crystal_models = candidate_orientation_matrices

