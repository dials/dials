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


class GaussianRSProfileModelExt(ProfileModelIface):
  ''' An extension class implementing a reciprocal space gaussian profile model. '''

  name = 'gaussian_rs'

  default=True

  @classmethod
  def phil(cls):
    from dials.algorithms.profile_model.gaussian_rs import phil_scope
    return phil_scope.as_str()

  def __init__(self):
    from dials.algorithms.profile_model.gaussian_rs import ProfileModelList
    self._model = ProfileModelList()

  def compute_bbox(self, experiments, reflections, **kwargs):
    return self._model.compute_bbox(experiments, reflections, **kwargs)

  def compute_partiality(self, experiments, reflections, **kwargs):
    return self._compute_partiality(experiments, reflections, **kwargs)

  def compute_mask(self, experiments, reflections, **kwargs):
    return self._model.compute_partiality(experiments, reflections, **kwargs)

  @classmethod
  def compute(cls, experiments, reflections, **kwargs):
    return ProfileModelList.compute(experiments, reflections, **kwargs)

  @classmethod
  def load(cls, params):
    return ProfileModelList.load(params)

  def dump(self):
    return self.dump()
