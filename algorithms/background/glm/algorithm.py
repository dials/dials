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
               model='constant3d',
               tuning_constant=1.345):
    '''
    Initialise the algorithm.

    :param experiments: The list of experiments
    :param model: The background model
    :param tuning_constant: The robust tuning constant

    '''
    from dials.algorithms.background.glm import Creator
    if model == 'constant2d':
      model = Creator.model.constant2d
    elif model == 'constant3d':
      model = Creator.model.constant3d
    elif model == 'loglinear2d':
      model = Creator.model.loglinear2d
    elif model == 'loglinear3d':
      model = Creator.model.loglinear3d
    else:
      raise RuntimeError('Unknown background model')
    self._create = Creator(
      model=model,
      tuning_constant=tuning_constant,
      max_iter=100)

  def compute_background(self, reflections):
    '''
    Compute the backgrond.

    :param reflections: The list of reflections

    '''
    from dials.array_family import flex

    # Do the background subtraction
    success = self._create(reflections['shoebox'])
    reflections['background.mean'] = flex.double(
      [sbox.background[0] for sbox in reflections['shoebox']])
    reflections.set_flags(success != True, reflections.flags.dont_integrate)
