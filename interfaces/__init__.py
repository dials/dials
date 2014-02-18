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

class CentroidIface(interface.Interface):
  ''' Interface for centroid calculation algorithms. '''

  name = 'centroid'

  @interface.abstractmethod
  def compute_centroid(self, reflections):
    pass


class BackgroundIface(interface.Interface):
  ''' Interface for background calculation algorithms. '''

  name = 'background'

  @interface.abstractmethod
  def compute_background(self, reflections):
    pass


class IntegrationIface(interface.Interface):
  ''' Interface for intensity calculation algorithms. '''
  name = 'integration'

  @interface.abstractmethod
  def compute_intensity(self, reflections):
    pass
