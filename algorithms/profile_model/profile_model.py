#!/usr/bin/env python
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


class ProfileModel(object):

  def __init__(self, experiment, reflections):
    self._sigma_b = 0
    self._sigma_m = 0

  def sigma_b(self):
    return self._sigma_b

  def sigma_m(self):
    return self._sigma_m
