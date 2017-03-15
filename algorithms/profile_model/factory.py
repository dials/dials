#
# factory.py
#
#  Copyright (C) 2013 Diamond Light Source
#
#  Author: James Parkhurst
#
#  This code is distributed under the BSD license, a copy of which is
#  included in the root directory of this package.
from __future__ import absolute_import, division

def generate_phil_scope():
  '''
  Generate the phil scope for profile model

  :return: The phil scope

  '''
  import dials.extensions
  from dials.interfaces import ProfileModelIface
  phil_scope = ProfileModelIface.phil_scope()
  return phil_scope

phil_scope = generate_phil_scope()

class ProfileModelFactory(object):
  '''
  Factory for creating profile models

  '''

  @classmethod
  def create(cls, params, experiments, reflections=None):
    '''
    Compute or load the profile model.

    :param params: The input phil parameters
    :param experiments: The experiment list
    :param reflections: The reflection table
    :return: The profile model

    '''
    from dials.interfaces import ProfileModelIface
    Extension = ProfileModelIface.extension(params.profile.algorithm)
    Algorithm = Extension().algorithm()
    create_profile_model = True
    if hasattr(params, 'create_profile_model'):
      create_profile_model = params.create_profile_model
    if (create_profile_model and
        reflections is not None and
        "shoebox" in reflections):
      for expr, indices in reflections.iterate_experiments_and_indices(experiments):
        expr.profile = Algorithm.create(
          params.profile,
          reflections.select(indices),
          expr.crystal,
          expr.beam,
          expr.detector,
          expr.goniometer,
          expr.scan)
    else:
      for expr in experiments:
        expr.profile = Algorithm.create(
          params.profile,
          None,
          expr.crystal,
          expr.beam,
          expr.detector,
          expr.goniometer,
          expr.scan,
          expr.profile)
    return experiments
