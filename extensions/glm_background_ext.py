#!/usr/bin/env python
#
# simple_background_ext.py
#
#  Copyright (C) 2013 Diamond Light Source
#
#  Author: James Parkhurst
#
#  This code is distributed under the BSD license, a copy of which is
#  included in the root directory of this package.
from __future__ import division

from dials.interfaces import BackgroundIface


class SimpleBackgroundExt(BackgroundIface):
  ''' An extension class implementing XDS background subtraction. '''

  name = 'glm'

  @classmethod
  def phil(cls):
    from libtbx.phil import parse
    phil = parse('''

      robust {
        tuning_constant = 1.345
          .type = float
          .help = "The tuning constant for robust estimation"
      }

      model {
        algorithm = constant2d *constant3d loglinear2d loglinear3d
          .type = choice
          .help = "The background model to fit"
      }

    ''')
    return phil

  def __init__(self, params, experiments):
    '''
    Initialise the algorithm.

    :param params: The input parameters
    :param experiments: The list of experiments

    '''
    from libtbx.phil import parse
    from dials.algorithms.background.glm import BackgroundAlgorithm

    # Create some default parameters
    if params is None:
      params = self.phil().fetch(parse('')).extract()
    else:
      params = params.integration.background.glm

    # Create the algorithm
    self._algorithm = BackgroundAlgorithm(
      experiments,
      tuning_constant=params.robust.tuning_constant,
      model=params.model.algorithm)

  def compute_background(self, reflections):
    '''
    Compute the background.

    :param reflections: The list of reflections

    '''
    self._algorithm.compute_background(reflections)
