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
  ''' Interface for threshold algorithms to be used in spot finding. '''

  scope = "spotfinder"
  name = 'threshold'

  def __init__(self, params, imageset):
    ''' Initialise the algorithm.

    :param params: The input phil parameters
    :param imageset: The imageset

    '''
    pass

  @interface.abstractmethod
  def compute_threshold(self, image, mask):
    ''' Threshold the image. The image may be either a flex.int or flex.double
    type. The thresholded image should be returned as a flex.bool where pixels
    labelled True are spot pixels and False are background.

    :param image: The image to threshold
    :param mask: The corresponding mask

    :returns: The thresholded image

    '''
    pass


class ProfileModelIface(interface.Interface):
  '''
  The interface definition for a profile model.

  '''

  name = 'profile'

  @interface.abstractmethod
  def algorithm(self):
    ''' Get the algorithm. '''
    pass


class CentroidIface(interface.Interface):
  ''' Interface for centroid algorithms. '''

  scope = "integration"
  name = 'centroid'

  def __init__(self, params, experiments):
    ''' Initialise the algorithm.

    :param params: The input phil parameters
    :param experiments: The experiment list

    '''
    pass

  @interface.abstractmethod
  def compute_centroid(self, reflections):
    ''' Compute the reflection centroids

    :param reflections: The list of reflections

    '''
    pass


class BackgroundIface(interface.Interface):
  ''' Interface for background algorithms. '''

  scope = "integration"
  name = 'background'

  def __init__(self, params, experiments):
    ''' Initialise the algorithm.

    :param params: The input phil parameters
    :param experiments: The experiment list

    '''
    pass

  @interface.abstractmethod
  def compute_background(self, reflections):
    ''' Compute the reflection background

    :param reflections: The list of reflections

    '''
    pass
