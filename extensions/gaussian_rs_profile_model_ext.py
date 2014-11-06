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

from dials.interfaces import ProfileModelCreatorIface


class GaussianRSProfileModelExt(ProfileModelCreatorIface):
  ''' An extension class implementing a reciprocal space gaussian profile model. '''

  name = 'gaussian_rs'

  default=True

  @classmethod
  def phil(cls):
    from dials.algorithms.profile_model.gaussian_rs import phil_scope
    return phil_scope

  def __init__(self):
    from dials.algorithms.profile_model.gaussian_rs import ProfileModelList
    self._model = ProfileModelList()

  @classmethod
  def create(cls, params, experiments, reflections):
    from dials.algorithms.profile_model.gaussian_rs import Factory
    return Factory.create(params, experiments, reflections)

