#!/usr/bin/env python
#
# gmodel_background_ext.py
#
#  Copyright (C) 2013 Diamond Light Source
#
#  Author: James Parkhurst
#
#  This code is distributed under the BSD license, a copy of which is
#  included in the root directory of this package.
from __future__ import division

from dials.interfaces import BackgroundIface


class GModelBackgroundExt(BackgroundIface):
  ''' An extension class implementing a robust GLM background algorithm. '''

  name = 'gmodel'

  @classmethod
  def phil(cls):
    from libtbx.phil import parse
    phil = parse('''

      robust {

        algorithm = False
          .type = bool
          .help = "Use the robust algorithm"

        tuning_constant = 1.345
          .type = float
          .help = "The tuning constant for robust estimation"
      }

      min_pixels = 10
        .type = int(value_min=1)
        .help = "The minimum number of pixels required"

      model = None
        .type = str
        .help = "The model filename"

    ''')
    return phil

  def __init__(self, params, experiments):
    '''
    Initialise the algorithm.

    :param params: The input parameters
    :param experiments: The list of experiments

    '''
    from libtbx.phil import parse
    from dials.algorithms.background.gmodel import BackgroundAlgorithm

    # Create some default parameters
    if params is None:
      params = self.phil().fetch(parse('')).extract()
    else:
      params = params.integration.background.gmodel

    # Create the algorithm
    self._algorithm = BackgroundAlgorithm(
      experiments,
      model           = params.model,
      robust          = params.robust.algorithm,
      tuning_constant = params.robust.tuning_constant,
      min_pixels      = params.min_pixels)

  def compute_background(self, reflections, image_volume=None):
    '''
    Compute the background.

    :param reflections: The list of reflections

    '''
    return self._algorithm.compute_background(
      reflections, image_volume=image_volume)
