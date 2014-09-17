#
# profile_model.py
#
#  Copyright (C) 2013 Diamond Light Source
#
#  Author: James Parkhurst
#
#  This code is distributed under the BSD license, a copy of which is
#  included in the root directory of this package.

from __future__ import division
from abc import ABCMeta, abstractmethod


class ProfileModelIface(object):
  ''' The interface definition for profile model algorithms. '''

  __metaclass__ = ABCMeta

  @abstractmethod
  def compute_bbox(self, experiment, reflections, **kwargs):
    ''' Given a single experiment and list of reflections, compute the bounding
    box of the reflections on the detector (and image frames).

    '''
    pass

  @abstractmethod
  def compute_partiality(self, experiment, reflections, **kwargs):
    ''' Given a single experiment and list of reflections, compute the
    partiality of the reflections

    '''
    pass

  @abstractmethod
  def compute_mask(self, experiment, reflections, **kwargs):
    ''' Given a single experiment and list of reflections, compute the
    foreground/background mask of the reflections.

    '''
    pass

  @abstractmethod
  def compute(cls, experiment, reflections, **kwargs):
    ''' Given an experiment and list of reflections, compute the profile model.

    This method should be a @classmethod and return an instance of the
    profile model class.

    '''
    pass


class ProfileModelListIface(object):
  ''' The interface definition for a list of profile models. '''

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
  def compute(cls, experiments, reflections, **kwargs):
    ''' Given a list of experiments and a list of reflections, compute the
    profile models.

    This method should be a @classmethod and return an instance of the
    profile model list class. '''
    pass

  @abstractmethod
  def load(cls, params):
    ''' Given a set of extracted phil parameters, load the profile model.

    This method should be a @classmethod and return an instance of the
    profile model list class. '''
    pass

  @abstractmethod
  def dump(self):
    ''' Dump and return the profile model to a phil scope object. '''
    pass
