#
# factory.py
#
#  Copyright (C) 2013 Diamond Light Source
#
#  Author: James Parkhurst
#
#  This code is distributed under the BSD license, a copy of which is
#  included in the root directory of this package.
from __future__ import division
from libtbx.phil import parse

def generate_phil_scope():
  import dials.extensions
  from dials.interfaces import ProfileModelIface
  phil_scope = ProfileModelIface.phil_scope()
  return phil_scope

phil_scope = generate_phil_scope()

class ProfileModelFactory(object):
  ''' Factory for creating profile models '''

  @classmethod
  def create(cls, params, experiments, reflections=None):
    ''' Compute or load the profile model. '''
    from dials.interfaces import ProfileModelIface
    Algorithm = ProfileModelIface.extension(params.profile.algorithm)
    return Algorithm.create(params, experiments, reflections)
