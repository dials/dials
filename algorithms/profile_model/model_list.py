#!/usr/bin/env python
#
# model_list.py
#
#  Copyright (C) 2013 Diamond Light Source
#
#  Author: James Parkhurst
#
#  This code is distributed under the BSD license, a copy of which is
#  included in the root directory of this package.
from __future__ import division

class ProfileModelList(object):
  '''
  A class to represent multiple profile models.

  '''

  def __init__(self, models=None):
    '''
    Initialise the model list.

    :param models: The profile models

    '''
    if models is None:
      self._models = []
    else:
      self._models = models

  def __getitem__(self, index):
    '''
    Get a profile model.

    :param index: The list index
    :return: A profile model or sublist

    '''
    if isinstance(index, slice):
      return ProfileModelList(self._models[index])
    return self._models[index]

  def __len__(self):
    '''
    :return: the number of models.

    '''
    return len(self._models)

  def __iter__(self):
    '''
    Iterate through the models.

    '''
    for i in range(len(self)):
      yield self[i]

  def append(self, model):
    '''
    Add another model.

    :param model: The profile model to append

    '''
    self._models.append(model)

  def predict_reflections(self, experiments, **kwargs):
    '''
    Predict the reflections.

    :param experiments: The experiment list
    :return: The predicted reflections

    '''
    from dials.array_family import flex
    result = flex.reflection_table()
    for index, (model, experiment) in enumerate(zip(self, experiments)):
      rlist = model.predict_reflections(experiment, **kwargs)
      rlist['id'] = flex.size_t(len(rlist), index)
      result.extend(rlist)
    return result

  def compute_bbox(self, experiments, reflections, **kwargs):
    '''
    Compute the bounding boxes.

    :param experiments: The experiment list
    :param reflections: The reflection table
    :return: The bounding boxes

    '''
    from dials.array_family import flex
    index_list = reflections.split_indices_by_experiment_id(len(experiments))
    assert(len(experiments) == len(self))
    assert(len(experiments) == len(index_list))
    reflections['bbox'] = flex.int6(len(reflections))
    for model, experiment, index in zip(self, experiments, index_list):
      reflections['bbox'].set_selected(
        index,
        model.compute_bbox(
          experiment,
          reflections.select(index),
          **kwargs))
    return reflections['bbox']

  def compute_partiality(self, experiments, reflections, **kwargs):
    '''
    Compute the partiality.

    :param experiments: The experiment list
    :param reflections: The reflection table
    :return: The partiality

    '''
    from dials.array_family import flex
    index_list = reflections.split_indices_by_experiment_id(len(experiments))
    assert(len(experiments) == len(self))
    assert(len(experiments) == len(index_list))
    reflections['partiality'] = flex.double(len(reflections))
    for model, experiment, index in zip(self, experiments, index_list):
      reflections['partiality'].set_selected(
        index,
        model.compute_partiality(
          experiment,
          reflections.select(index),
          **kwargs))
    return reflections['partiality']

  def compute_mask(self, experiments, reflections, **kwargs):
    '''
    Compute the shoebox mask.

    :param experiments: The experiment list
    :param reflections: The reflection table

    '''
    from dials.array_family import flex
    index_list = reflections.split_indices_by_experiment_id(len(experiments))
    assert(len(experiments) == len(self))
    assert(len(experiments) == len(index_list))
    for model, experiment, index in zip(self, experiments, index_list):
      model.compute_mask(
        experiment,
        reflections.select(index),
        **kwargs)

  def dump(self):
    '''
    Dump the profile model to phil parameters.

    :return: The phil parameters

    '''
    from dials.algorithms.profile_model import factory
    return factory.phil_scope.fetch(sources=[model.dump() for model in self])
