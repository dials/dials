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

from __future__ import division
import math

from scitbx import matrix
from scitbx.array_family import flex
from dials.algorithms.indexing.indexer2 import \
     indexer_base, optimise_basis_vectors
from dials.algorithms.indexing.indexer2 import \
     is_approximate_integer_multiple



class indexer_real_space_grid_search(indexer_base):

  def __init__(self, reflections, sweep, params):
    super(indexer_real_space_grid_search, self).__init__(
      reflections, sweep, params)

  def find_lattices(self):
    return self.find_lattices_real_space_grid_search()

  def find_lattices_real_space_grid_search(self):
    self.d_min = self.params.refinement_protocol.d_min_start # XXX this should be another parameter
    self.real_space_grid_search()
    crystal_models = self.candidate_crystal_models
    return crystal_models

  def real_space_grid_search(self):
    d_min = self.params.refinement_protocol.d_min_start

    reciprocal_space_points = self.reciprocal_space_points.select(
      (self.reflections.crystal() == -1) &
      (1/self.reciprocal_space_points.norms() > d_min))
    #reciprocal_space_points = reciprocal_space_points.select(
      #1/reciprocal_space_points.norms() > self.params.refinement_protocol.d_min_start)

    print "indexing from %i reflections" %len(reciprocal_space_points)

    def compute_functional(vector):
      two_pi_S_dot_v = 2 * math.pi * reciprocal_space_points.dot(vector)
      return flex.sum(flex.cos(two_pi_S_dot_v))

    from rstbx.array_family import flex
    from rstbx.dps_core import SimpleSamplerTool
    assert self.target_symmetry_primitive is not None
    assert self.target_symmetry_primitive.unit_cell() is not None
    SST = SimpleSamplerTool(0.020)
    SST.construct_hemisphere_grid(SST.incr)
    cell_dimensions = self.target_symmetry_primitive.unit_cell().parameters()[:3]
    unique_cell_dimensions = set(cell_dimensions)
    print "Number of search vectors: %i" %(len(SST.angles) * len(unique_cell_dimensions))
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

    #for i in range(30):
      #v = matrix.col(vectors[i])
      #print v.elems, v.length(), function_values[i]

    basis_vectors = [v.elems for v in unique_vectors]
    self.candidate_basis_vectors = basis_vectors

    if self.params.optimise_initial_basis_vectors:
      optimised_basis_vectors = optimise_basis_vectors(
        reciprocal_space_points, basis_vectors, n_multiples=1)
      optimised_function_values = flex.double([
        compute_functional(v) for v in optimised_basis_vectors])

      perm = flex.sort_permutation(optimised_function_values, reverse=True)
      optimised_basis_vectors = optimised_basis_vectors.select(perm)
      optimised_function_values = optimised_function_values.select(perm)

      unique_vectors = [matrix.col(v) for v in optimised_basis_vectors]

    print "Number of unique vectors: %i" %len(unique_vectors)

    for i in range(len(unique_vectors)):
      print compute_functional(unique_vectors[i].elems), unique_vectors[i].length(), unique_vectors[i].elems
      print

    crystal_models = []
    while True:
      self.candidate_basis_vectors = unique_vectors
      if self.params.debug:
        self.debug_show_candidate_basis_vectors()
      candidate_orientation_matrices \
        = self.find_candidate_orientation_matrices(
          unique_vectors, return_first=False)
      if len(candidate_orientation_matrices) == 0: break
      n_indexed = flex.int()
      from dials.algorithms.indexing import index_reflections
      for cm in candidate_orientation_matrices:
        refl = self.reflections.deep_copy().select(
          (self.reflections.crystal() == -1) &
          (1/self.reciprocal_space_points.norms() > d_min))
        refl.set_crystal(flex.int(refl.size(), -1))
        index_reflections(refl, reciprocal_space_points,
                          [cm], self.d_min, tolerance=0.25,
                          verbosity=0)
        n_indexed.append((refl.crystal() > -1).count(True))
      perm = flex.sort_permutation(n_indexed, reverse=True)
      print list(perm)
      print list(n_indexed.select(perm))
      print candidate_orientation_matrices[perm[0]]

      if n_indexed[perm[0]] < 50:
        break

      new_unique_vectors = []
      cm = candidate_orientation_matrices[perm[0]]
      crystal_models.append(cm)
      a, b, c = cm.get_real_space_vectors()
      for v in unique_vectors:
        if not (is_approximate_integer_multiple(v, a) or
                is_approximate_integer_multiple(v, b) or
                is_approximate_integer_multiple(v, c)):
          new_unique_vectors.append(v)
      assert len(new_unique_vectors) == len(unique_vectors) - 3
      if len(new_unique_vectors) < 3: break
      unique_vectors = new_unique_vectors
      if len(crystal_models) == 1:
        break

    #assert len(crystal_models) > 0

    candidate_orientation_matrices = crystal_models

    for i in range(len(candidate_orientation_matrices)):
      if self.target_symmetry_primitive is not None:
        #print "symmetrizing model"
        #self.target_symmetry_primitive.show_summary()
        symmetrized_model = self.apply_symmetry(
          candidate_orientation_matrices[i], self.target_symmetry_primitive)
        candidate_orientation_matrices[i] = symmetrized_model

    self.candidate_crystal_models = candidate_orientation_matrices

