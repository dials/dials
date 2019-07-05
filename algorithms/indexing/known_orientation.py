#!/usr/bin/env python
# -*- mode: python; coding: utf-8; indent-tabs-mode: nil; python-indent: 2 -*-
#
# dials.algorithms.indexing.known_orientation.py
#
#  Copyright (C) 2014 Diamond Light Source
#
#  Author: Richard Gildea
#
#  This code is distributed under the BSD license, a copy of which is
#  included in the root directory of this package.

from __future__ import absolute_import, division, print_function

from dials.algorithms.indexing.indexer import Indexer
from dxtbx.model.experiment_list import Experiment, ExperimentList


class IndexerKnownOrientation(Indexer):
    def __init__(self, reflections, experiments, params, known_orientations):
        self.known_orientations = known_orientations
        super(IndexerKnownOrientation, self).__init__(reflections, experiments, params)

    def find_lattices(self):
        experiments = ExperimentList()
        for cm in self.known_orientations:
            # indexer expects crystals to be in primitive setting
            space_group = cm.get_space_group()
            cb_op_to_primitive = (
                space_group.info().change_of_basis_op_to_primitive_setting()
            )
            cm = cm.change_basis(cb_op_to_primitive)
            for expt in self.experiments:
                experiments.append(
                    Experiment(
                        imageset=expt.imageset,
                        beam=expt.beam,
                        detector=expt.detector,
                        goniometer=expt.goniometer,
                        scan=expt.scan,
                        crystal=cm,
                    )
                )
        return experiments
