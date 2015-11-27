#!/usr/bin/env python
#
# const_d_background_ext.py
#
#  Copyright (C) 2013 Diamond Light Source
#
#  Author: James Parkhurst
#
#  This code is distributed under the BSD license, a copy of which is
#  included in the root directory of this package.
from __future__ import division

from dials.interfaces import BackgroundIface


class ConstDBackgroundExt(BackgroundIface):
  ''' An extension class implementing a crude background algorithm. '''

  name = 'const_d'

  @classmethod
  def phil(cls):
    from libtbx.phil import parse
    phil = parse('''

    ''')
    return phil

  def __init__(self, params, experiments):
    '''
    Initialise the algorithm.

    :param params: The input parameters
    :param experiments: The list of experiments

    '''
    from libtbx.phil import parse
    from dials.algorithms.background.const_d import BackgroundAlgorithm

    # Create the algorithm
    self._algorithm = BackgroundAlgorithm(experiments)

  def compute_background(self, reflections):
    '''
    Compute the background.

    :param reflections: The list of reflections

    '''
    return self._algorithm.compute_background(reflections)
