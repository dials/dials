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
from dxtbx.model.profile import Profile as ProfileBase
from abc import abstractmethod


class ProfileModelIface(ProfileBase):
  '''
  The interface definition for a list of profile models.

  '''

  @abstractmethod
  def predict_reflections(self, experiment, **kwargs):
    '''
    Given an experiment, predict the reflections.

    :param experiment: The experiment

    '''
    pass

  @abstractmethod
  def compute_partiality(self, experiment, reflections, **kwargs):
    '''
    Given an experiment and list of reflections, compute the partiality of the
    reflections

    :param experiment: The experiment
    :param reflections: The reflection table

    '''
    pass

  @abstractmethod
  def compute_bbox(self, experiment, reflections, **kwargs):
    ''' Given an experiment and list of reflections, compute the
    bounding box of the reflections on the detector (and image frames).

    '''
    pass

  @abstractmethod
  def compute_mask(self, experiment, reflections, **kwargs):
    '''
    Given an experiment and list of reflections, compute the
    foreground/background mask of the reflections.

    :param experiment: The experiment
    :param reflections: The reflection table

    '''
    pass

