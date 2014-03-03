#!/usr/bin/env python
#
# __init__.py
#
#  Copyright (C) 2013 Diamond Light Source
#
#  Author: James Parkhurst
#
#  This code is distributed under the BSD license, a copy of which is
#  included in the root directory of this package.

from __future__ import division
from dials.framework import interface


class SpotFinderThresholdIface(interface.Interface):
  ''' Configure the spot finder threshold algorithms. '''

  name = 'threshold'

  def __init__(self, params, imageset):
    pass

  @interface.abstractmethod
  def compute_threshold(self, image):
    pass


class CentroidIface(interface.Interface):
  ''' Configure the centroid calculation algorithms. '''

  name = 'centroid'

  def __init__(self, params, experiment):
    pass

  @interface.abstractmethod
  def compute_centroid(self, reflections):
    pass


class BackgroundIface(interface.Interface):
  ''' Configure the background calculation algorithms. '''

  name = 'background'

  def __init__(self, params, experiment):
    pass

  @interface.abstractmethod
  def compute_background(self, reflections):
    pass


class IntegrationIface(interface.Interface):
  ''' Configure the intensity calculation algorithms. '''

  name = 'integration'

  def __init__(self, params, experiment):
    pass

  @interface.abstractmethod
  def compute_intensity(self, reflections):
    pass
