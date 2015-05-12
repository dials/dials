#!/usr/bin/env python
#
# general.py
#
#  Copyright (C) 2013 Diamond Light Source
#
#  Author: James Parkhurst
#
#  This code is distributed under the BSD license, a copy of which is
#  included in the root directory of this package.

from __future__ import division

class BackgroundAlgorithm(object):
  ''' Class to do background subtraction. '''

  def __init__(self, experiments,
               tuning_constant=1.345):
    '''
    Initialise the algorithm.

    :param experiments: The list of experiments
    :param tuning_constant: The robust tuning constant

    '''
    from dials.algorithms.background.glm import Modeller
    self._model = Modeller(
      tuning_constant=tuning_constant,
      max_iter=100)

  def compute_background(self, reflections):
    '''
    Compute the backgrond.

    :param reflections: The list of reflections

    '''
    from dials.array_family import flex

    # Do the background subtraction
    success = self._model(reflections['shoebox'])
    reflections['background.mean'] = reflections['shoebox'].background[0]
    reflections.set_flags(success != True, reflections.flags.dont_integrate)
