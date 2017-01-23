#
# profile.py
#
#  Copyright (C) 2015 Diamond Light Source
#
#  Author: James Parkhurst
#
#  This code is distributed under the BSD license, a copy of which is
#  included in the root directory of this package.

from __future__ import absolute_import, division

from dxtbx.model.profile import ProfileModelBaseIface
import abc

class ProfileModelIface(ProfileModelBaseIface):
  '''
  The abc definition for a profile model.

  '''

  @classmethod
  def create(cls,
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

  @abc.abstractmethod
  def predict_reflections(self,
                          imageset,
                          crystal,
                          beam,
                          detector,
                          goniometer=None,
                          scan=None,
                          **kwargs):
    '''
    Given an experiment, predict the reflections.

    :param imageset: The imageset
    :param crystal: The crystal model
    :param beam: The beam model
    :param detector: The detector model
    :param goniometer: The goniometer model
    :param scan: The scan model

    '''
    pass

  @abc.abstractmethod
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

  @abc.abstractmethod
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

  @abc.abstractmethod
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

  def fitting_class(self):
    '''
    Get the profile fitting algorithm associated with this profile model

    :return: The profile fitting class

    '''
    return None
