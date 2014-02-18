#!/usr/bin/env python
#
# dials.algorithms.centroid.centroid_factory.py
#
#  Copyright (C) 2013 Diamond Light Source
#
#  Author: James Parkhurst
#
#  This code is distributed under the BSD license, a copy of which is
#  included in the root directory of this package.

from __future__ import division


class CentroidFactory(object):
  ''' Factory class to create integrators '''

  @staticmethod
  def from_parameters(params):
    ''' Given a set of parameters, configure the centroid calculator

    Params:
        params The input parameters

    Returns:
        The centroid calculator instance

    '''
    from dials.algorithms.centroid.centroider import Centroider
    return Centroider()

