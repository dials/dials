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

class Algorithm(object):
  '''
  A python class to wrap the centroid algorithm

  '''

  def __init__(self, experiments):
    '''
    Initialize the centroider

    :param experiments: The experiment list

    '''
    from dials.algorithms.centroid.simple import Centroider

    # Create the centroider
    self.centroider = Centroider()

    # Add all the experiments
    for exp in experiments:
      if exp.scan is not None:
        self.centroider.add(
          exp.detector,
          exp.scan)
      else:
        self.centroider.add(
          exp.detector)

  def __call__(self, reflections, image_volume=None):
    '''
    Do the centroiding

    :param reflections: The reflection list

    '''
    if image_volume is None:
      return self.centroider(reflections)
    return self.centroider(reflections, image_volume)
