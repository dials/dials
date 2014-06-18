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

from dials.algorithms.indexing.indexer2 import indexer_base
from dxtbx.model.experiment.experiment_list import Experiment, ExperimentList

class indexer_fft1d(indexer_base):

  def __init__(self, reflections, sweep, params):
    super(indexer_fft1d, self).__init__(reflections, sweep, params)

  def find_lattices(self):
    self.d_min = self.params.refinement_protocol.d_min_start

    from dials.algorithms.indexing.indexer import candidate_basis_vectors_fft1d

    from rstbx.phil.phil_preferences import indexing_api_defs
    import iotbx.phil
    hardcoded_phil = iotbx.phil.parse(
      input_string=indexing_api_defs).extract()

    reflections = self.reflections.select(
      (self.reflections['id'] == -1) &
      (1/self.reciprocal_space_points.norms() > self.d_min))
    solutions = candidate_basis_vectors_fft1d(
      reflections, self.detector, self.beam,
      self.goniometer, self.scan, hardcoded_phil,
      max_cell=self.params.max_cell)
    self.candidate_basis_vectors = solutions[0]
    if self.params.debug:
      self.debug_show_candidate_basis_vectors()
    self.candidate_crystal_models = self.find_candidate_orientation_matrices(
      self.candidate_basis_vectors, return_first=True, apply_symmetry=False)
    crystal_models = self.candidate_crystal_models[:1]
    experiments = ExperimentList()
    for cm in crystal_models:
      experiments.append(Experiment(beam=self.beam,
                                    detector=self.detector,
                                    goniometer=self.goniometer,
                                    scan=self.scan,
                                    crystal=cm))
    return experiments
