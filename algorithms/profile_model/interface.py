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
  def compute_bbox(self, experiments, reflections, **kwargs):
    ''' Given a list of experiments and list of reflections, compute the
    bounding box of the reflections on the detector (and image frames).

    '''
    pass

  @abstractmethod
  def compute_partiality(self, experiments, reflections, **kwargs):
    ''' Given a list of experiments and list of reflections, compute the
    partiality of the reflections

    '''
    pass

  @abstractmethod
  def compute_mask(self, experiments, reflections, **kwargs):
    ''' Given a list of experiments and list of reflections, compute the
    foreground/background mask of the reflections.

    '''
    pass

  @abstractmethod
  def dump(self):
    ''' Dump and return the profile model to a phil scope object. '''
    pass

  @abstractmethod
  def __len__(self):
    ''' The number of models (should equal the number of experiments. '''
    pass

  @abstractmethod
  def __iter__(self):
    ''' Iterate through models '''
    pass

  @abstractmethod
  def __getitem__(self, index):
    ''' Get a model. '''
    pass
