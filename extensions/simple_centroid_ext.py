#!/usr/bin/env python
#
# simple_centroid_ext.py
#
#  Copyright (C) 2013 Diamond Light Source
#
#  Author: James Parkhurst
#
#  This code is distributed under the BSD license, a copy of which is
#  included in the root directory of this package.
from __future__ import absolute_import, division

from dials.interfaces import CentroidIface


class SimpleCentroidExt(CentroidIface):
  ''' An extension class implementing a simple centroid algorithm. '''

  name = 'simple'

  default = True

  def __init__(self, params, experiments):
    ''' Initialise the algorithm.

    :param params: The input phil parameters
    :param experiments: The experiment list

    '''
    self.experiments = experiments

  def compute_centroid(self, reflections, image_volume=None):
    '''
    Compute the centroid.

    :param reflections: The list of reflections

    '''
    from dials.algorithms.centroid.simple.algorithm import Algorithm
    algorithm = Algorithm(self.experiments)
    return algorithm(reflections, image_volume=image_volume)
