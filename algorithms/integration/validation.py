#
# validation.py
#
#  Copyright (C) 2013 Diamond Light Source
#
#  Author: James Parkhurst
#
#  This code is distributed under the BSD license, a copy of which is
#  included in the root directory of this package.

from __future__ import absolute_import, division

class ValidatedMultiExpProfileModeller(object):
  '''
  A class to wrap profile modeller for validation

  '''

  def __init__(self):
    '''
    Init the list of modellers

    '''
    self.modellers = []
    self.finalized_modeller = None

  def add(self, modeller):
    '''
    Add a MultiExpProfileModeller

    '''
    self.modellers.append(modeller)

  def __getitem__(self, index):
    '''
    Get a modeller

    '''
    return self.modellers[index]

  def model(self, reflections):
    '''
    Do the modelling for all modellers

    '''
    from dials.array_family import flex
    if 'profile.index' not in reflections:
      assert(len(self.modellers) == 1)
      self.modellers[0].model(reflections)
    else:
      for i, modeller in enumerate(self.modellers):
        mask = reflections['profile.index'] != i
        indices = flex.size_t(range(len(mask))).select(mask)
        if len(indices) > 0:
          subsample = reflections.select(indices)
          modeller.model(subsample)
          reflections.set_selected(indices, subsample)

  def validate(self, reflections):
    '''
    Do the validation.

    '''
    from dials.array_family import flex
    results = []
    for i, modeller in enumerate(self.modellers):
      mask = reflections['profile.index'] != i
      indices = flex.size_t(range(len(mask))).select(mask)
      if len(indices) > 0:
        subsample = reflections.select(indices)
        modeller.validate(subsample)
        reflections.set_selected(indices, subsample)
        corr = subsample['profile.correlation']
        mean_corr = flex.mean(corr)
      else:
        mean_corr = None
      results.append(mean_corr)
    return results

  def accumulate(self, other):
    '''
    Accumulate the modellers

    '''
    assert(len(self) == len(other))
    for ms, mo in zip(self, other):
      ms.accumulate(mo)

  def finalize(self):
    '''
    Finalize the model

    '''
    assert not self.finalized()
    for m in self:
      if self.finalized_modeller is None:
        self.finalized_modeller = m.copy()
      else:
        self.finalized_modeller.accumulate(m)
      m.finalize()
    self.finalized_modeller.finalize()

  def finalized(self):
    '''
    Check if the model has been finalized.

    '''
    return self.finalized_modeller is not None

  def finalized_model(self):
    '''
    Get the finalized model

    '''
    assert self.finalized
    return self.finalized_modeller

  def __iter__(self):
    '''
    Iterate through the modellers

    '''
    for m in self.modellers:
      yield m

  def __len__(self):
    '''
    Return the number of modellers

    '''
    return len(self.modellers)
