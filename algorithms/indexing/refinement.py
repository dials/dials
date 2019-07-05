#!/usr/bin/env python
# -*- mode: python; coding: utf-8; indent-tabs-mode: nil; python-indent: 2 -*-
# dials.algorithms.indexing.refinement.py
#
#  Copyright (C) 2014 Diamond Light Source
#
#  Author: Richard Gildea
#
#  This code is distributed under the BSD license, a copy of which is
#  included in the root directory of this package.

from __future__ import absolute_import, division, print_function


import logging

logger = logging.getLogger(__name__)


def refine(params, reflections, experiments, verbosity=0):
    if params.refinement.parameterisation.scan_varying:
        logger.warning(
            "scan_varying=True not supported in indexing: setting scan_varying=False"
        )
        params.refinement.parameterisation.scan_varying = False

    from dials.algorithms.refinement import RefinerFactory

    refiner = RefinerFactory.from_parameters_data_experiments(
        params, reflections, experiments, verbosity=verbosity
    )

    outliers = None
    refined = refiner.run()
    return refiner, refined, outliers
