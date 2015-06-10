#!/usr/bin/env python
#
# __init__.py
#
#  Copyright (C) 2015 Diamond Light Source
#
#  Author: James Parkhurst
#
#  This code is distributed under the BSD license, a copy of which is
#  included in the root directory of this package.

from __future__ import division
from dxtbx.model.profile import ProfileModelBaseIface
from abc import abstractmethod


class ProfileModelIface(ProfileBaseIface):
  '''
  The interface definition for a list of profile models.

  '''

  @classmethod
  def create(Class, params, experiment, reflections):
    '''
    Create the profile model from data.

    :param params: The phil parameters
    :param experiment: The experiment
    :param reflections: The reflections
    :return: An instance of the profile model

    '''
    return None

  @abstractmethod
  def predict_reflections(self, experiment):
    '''
    Given an experiment, predict the reflections.

    :param experiment: The experiment

    '''
    pass

  @abstractmethod
  def compute_partiality(self, experiment, reflections):
    '''
    Given an experiment and list of reflections, compute the partiality of the
    reflections

    :param experiment: The experiment
    :param reflections: The reflection table

    '''
    pass

  @abstractmethod
  def compute_bbox(self, experiment, reflections):
    ''' Given an experiment and list of reflections, compute the
    bounding box of the reflections on the detector (and image frames).

    :param experiment: The experiment
    :param reflections: The reflection table

    '''
    pass

  @abstractmethod
  def compute_mask(self, experiment, reflections):
    '''
    Given an experiment and list of reflections, compute the
    foreground/background mask of the reflections.

    :param experiment: The experiment
    :param reflections: The reflection table

    '''
    pass

  @classmethod
  def profile_fitting_class(Class):
    '''
    Get the profile fitting algorithm associated with this profile model

    :return: The profile fitting class

    '''
    return None
