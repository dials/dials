#!/usr/bin/env python
#
# gaussian_rs_profile_model_ext.py
#
#  Copyright (C) 2013 Diamond Light Source
#
#  Author: James Parkhurst
#
#  This code is distributed under the BSD license, a copy of which is
#  included in the root directory of this package.
from __future__ import division

from dials.interfaces import ProfileModelIface
from dials.algorithms.profile_model.gaussian_rs import Model


class GaussianRSProfileModelExt(ProfileModelIface):
  ''' An extension class implementing a reciprocal space gaussian profile model. '''

  name = 'gaussian_rs'

  default=True

  @classmethod
  def phil(cls):
    from dials.algorithms.profile_model.gaussian_rs import phil_scope
    return phil_scope

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
    return Model.create(
      params,
      reflections,
      crystal,
      beam,
      detector,
      goniometer,
      scan)

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
    return Model.predict_reflections(
      crystal,
      beam,
      detector,
      goniometer,
      scan,
      **kwargs)

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
    return Model.compute_partiality(
      reflections,
      crystal,
      beam,
      detector,
      goniometer,
      scan,
      **kwargs)

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
    return Model.compute_bbox(
      reflections,
      crystal,
      beam,
      detector,
      goniometer,
      scan,
      **kwargs)

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
    return Model.compute_mask(
      reflections,
      crystal,
      beam,
      detector,
      goniometer,
      scan,
      **kwargs)

  @classmethod
  def profile_fitting_class(Class):
    '''
    Get the profile fitting algorithm associated with this profile model

    :return: The profile fitting class

    '''
    return None
