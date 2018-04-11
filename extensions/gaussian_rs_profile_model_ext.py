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
from __future__ import absolute_import, division

class GaussianRSProfileModelExt(object):
  ''' An extension class implementing a reciprocal space gaussian profile model. '''

  name = 'gaussian_rs'

  default=True

  @classmethod
  def phil(cls):
    from dials.algorithms.profile_model.gaussian_rs import phil_scope
    return phil_scope

  @classmethod
  def algorithm(cls):
    from dials.algorithms.profile_model.gaussian_rs import Model
    return Model

  @classmethod
  def from_dict(cls, d):
    return cls.algorithm().from_dict(d)
