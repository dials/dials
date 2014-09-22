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
  ''' The interface definition for a list of profile models. '''

  name = 'profile'

  @interface.abstractmethod
  def compute_bbox(self, experiments, reflections, **kwargs):
    ''' Given a list of experiments and list of reflections, compute the
    bounding box of the reflections on the detector (and image frames).

    '''
    pass

  @interface.abstractmethod
  def compute_partiality(self, experiments, reflections, **kwargs):
    ''' Given a list of experiments and list of reflections, compute the
    partiality of the reflections

    '''
    pass

  @interface.abstractmethod
  def compute_mask(self, experiments, reflections, **kwargs):
    ''' Given a list of experiments and list of reflections, compute the
    foreground/background mask of the reflections.

    '''
    pass

  @interface.abstractmethod
  def dump(self):
    ''' Dump and return the profile model to a phil scope object. '''
    pass

  @interface.abstractmethod
  def create(cls, params, experiments, reflections):
    ''' Given a list of experiments and a list of reflections, compute the
    profile models.

    This method should be a @classmethod and return an instance of the
    profile model list class. '''
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


class IntensityIface(interface.Interface):
  ''' Interface for intensity algorithms. '''

  scope = "integration"
  name = 'intensity'

  def __init__(self, params, experiments, profile_model):
    ''' Initialise the algorithm.

    :param params: The input phil parameters
    :param experiments: The experiment list

    '''
    pass

  @interface.abstractmethod
  def type(self, params, experiments):
    ''' Return the type of integrator.

    This function should be implemented as a classmethod. The function provides
    information about how the shoeboxes should be read from file. The return
    value should be one of the following values:

     - 3d - shoeboxes will be read in 3D.
     - flat3d - shoeboxes will be read in 3D and summed across all frames.
     - 2d - shoeboxes will be read as single frame partials.
     - single2d - shoeboxes will be read for a single frame at a time.
     - stills - shoeboxes will be read for a single frame at a time for stills

    :param params: The input phil parameters
    :param experiments: The experiment list

    :returns: The type of integrator needed by this algorithm.

    '''
    pass

  @interface.abstractmethod
  def compute_intensity(self, reflections):
    ''' Compute the reflection intensity

    :param reflections: The list of reflections

    '''
    pass
