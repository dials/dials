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

phil_scope = parse('''
  profile {

    algorithm = *gaussian_rs
      .type = choice
      .help = "The profile model to use"

    include scope dials.algorithms.profile_model.gaussian_rs.phil_scope

  }
''', process_includes=True)


class ProfileModelFactory(object):
  ''' Factory for creating profile models '''

  @classmethod
  def create(cls, params, experiments, reflections, **kwargs):
    ''' Compute or load the profile model. '''

    # Select the algorithm and get the profile model
    if params.profile.algorithm == 'gaussian_rs':
      if len(params.profile.gaussian_rs) > 0:
        assert(len(params.profile.gaussian_rs) == len(experiments))
        model = ProfileModelFactory.load(params)
      else:
        assert(reference is not None)
        model = ProfileModelFactory.compute(params, experiment, reflections)
    else:
      raise RuntimeError(
        'Unknown profile algorithm %s' % (
          params.profile.algorithm))

    # Return the model
    return model

  @classmethod
  def compute(cls, params, experiments, reflections, **kwargs):
    ''' Compute the profile model. '''
    from dials.algorithms.profile_model import gaussian_rs

    # Select the algorithm and compute the profile model
    if params.profile.algorithm == 'gaussian_rs':
      model = gaussian_rs.ProfileModelList.compute(
        experiments,
        reflections)
    else:
      raise RuntimeError(
        'Unknown profile algorithm %s' % (
          params.profile.algorithm))

    # Return the profile model
    return model

  @classmethod
  def load(cls, params):
    ''' Load the profile model. '''
    from dials.algorithms.profile_model import gaussian_rs

    # Select the algorithm and compute the profile model
    if params.profile.algorithm == 'gaussian_rs':
      model = gaussian_rs.ProfileModelList.load(
        params.profile)
    else:
      raise RuntimeError(
        'Unknown profile algorithm %s' % (
          params.profile.algorithm))

    # Return the profile model
    return model
