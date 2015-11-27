#!/usr/bin/env python
#
# algorithm.py
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

  def __init__(self, experiments):
    '''
    Initialise the algorithm.

    :param experiments: The list of experiments

    '''
    from dials.algorithms.background.const_d import Creator
    assert len(experiments) == 1
    self._create = Creator(
      experiments[0].beam,
      experiments[0].detector)

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
    return success
