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
from dxtbx.model.profile import ProfileModelBaseIface


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


class ProfileModelIface(interface.Interface, ProfileBaseIface):
  '''
  The interface definition for a profile model.

  '''

  name = 'profile'

  @classmethod
  def create(Class,
             params,
             reflections,
             crystal,
             beam,
             detector,
             goniometer=None,
             scan=None):
    '''
    Create the profile model from data.

    :param params: The phil parameters
    :param reflections: The reflections
    :param crystal: The crystal model
    :param beam: The beam model
    :param detector: The detector model
    :param goniometer: The goniometer model
    :param scan: The scan model
    :return: An instance of the profile model

    '''
    return None

  @interface.abstractmethod
  def predict_reflections(self,
                          crystal,
                          beam,
                          detector,
                          goniometer=None,
                          scan=None,
                          **kwargs):
    '''
    Given an experiment, predict the reflections.

    :param crystal: The crystal model
    :param beam: The beam model
    :param detector: The detector model
    :param goniometer: The goniometer model
    :param scan: The scan model

    '''
    pass

  @interface.abstractmethod
  def compute_partiality(self,
                         reflections,
                         crystal,
                         beam,
                         detector,
                         goniometer=None,
                         scan=None,
                         **kwargs):
    '''
    Given an experiment and list of reflections, compute the partiality of the
    reflections

    :param reflections: The reflection table
    :param crystal: The crystal model
    :param beam: The beam model
    :param detector: The detector model
    :param goniometer: The goniometer model
    :param scan: The scan model

    '''
    pass

  @interface.abstractmethod
  def compute_bbox(self,
                   reflections,
                   crystal,
                   beam,
                   detector,
                   goniometer=None,
                   scan=None,
                   **kwargs):
    ''' Given an experiment and list of reflections, compute the
    bounding box of the reflections on the detector (and image frames).

    :param reflections: The reflection table
    :param crystal: The crystal model
    :param beam: The beam model
    :param detector: The detector model
    :param goniometer: The goniometer model
    :param scan: The scan model

    '''
    pass

  @interface.abstractmethod
  def compute_mask(self,
                   reflections,
                   crystal,
                   beam,
                   detector,
                   goniometer=None,
                   scan=None,
                   **kwargs):
    '''
    Given an experiment and list of reflections, compute the
    foreground/background mask of the reflections.

    :param reflections: The reflection table
    :param crystal: The crystal model
    :param beam: The beam model
    :param detector: The detector model
    :param goniometer: The goniometer model
    :param scan: The scan model

    '''
    pass

  @classmethod
  def profile_fitting_class(Class):
    '''
    Get the profile fitting algorithm associated with this profile model

    :return: The profile fitting class

    '''
    return None


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
