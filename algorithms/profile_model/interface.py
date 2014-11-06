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
from abc import abstractmethod, ABCMeta


class ProfileModelIface(object):
  ''' The interface definition for a list of profile models. '''

  __metaclass__ = ABCMeta

  @abstractmethod
  def predict_reflections(self, experiment, **kwargs):
    ''' Given an experiment, predict the reflections. '''
    pass

  @abstractmethod
  def compute_partiality(self, experiment, reflections, **kwargs):
    ''' Given an experiment and list of reflections, compute the
    partiality of the reflections

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
    ''' Given an experiment and list of reflections, compute the
    foreground/background mask of the reflections.

    '''
    pass

  @abstractmethod
  def dump(self):
    ''' Dump and return the profile model to a phil scope object. '''
    pass

