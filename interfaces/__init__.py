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

  scope = "spotfinder"
  name = 'threshold'

  def __init__(self, params, imageset):
    pass

  @interface.abstractmethod
  def compute_threshold(self, image, mask):
    pass


class CentroidIface(interface.Interface):
  ''' Configure the centroid calculation algorithms. '''

  scope = "integration"
  name = 'centroid'

  def __init__(self, params, experiments):
    pass

  @interface.abstractmethod
  def compute_centroid(self, reflections):
    pass


class BackgroundIface(interface.Interface):
  ''' Configure the background calculation algorithms. '''

  scope = "integration"
  name = 'background'

  def __init__(self, params, experiments):
    pass

  @interface.abstractmethod
  def compute_background(self, reflections):
    pass


class IntensityIface(interface.Interface):
  ''' Configure the intensity calculation algorithms. '''

  scope = "integration"
  name = 'intensity'

  def __init__(self, params, experiments, profile_model):
    pass

  @interface.abstractmethod
  def compute_intensity(self, reflections):
    pass


class Integration3DMixin(object):
  ''' A mixin to identify 3D integration algorithms. '''
  pass


class IntegrationFlat2DMixin(object):
  ''' A mixin to identify flat 2D integration algorithms. '''
  pass


class Integration2DMixin(object):
  ''' A mixin to identify 2D integration algorithms. '''
  pass


class IntegrationStillsMixin(object):
  ''' A mixin to identify 2D stills integration algorithms. '''
  pass
