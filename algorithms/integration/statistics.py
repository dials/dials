#
# statistics.py
#
#  Copyright (C) 2013 Diamond Light Source
#
#  Author: James Parkhurst
#
#  This code is distributed under the BSD license, a copy of which is
#  included in the root directory of this package.
from __future__ import division
from dials.array_family import flex
from dials.array_family.flex import Binner


class ModellerSummary(object):
  ''' A class to present a summary of integration results. '''

  def __init__(self,
               index,
               data,
               experiment,
               model
               ):
    ''' Initialise. '''
    ntot = len(model)
    ngood = 0
    nbad = 0
    for i in range(ntot):
      try:
        data = model.data(i)
        ngood += 1
      except Exception:
        nbad += 1
    self.text  = "Summary of profile modelling:\n"
    self.text += " Tried to model %d profiles\n" % ntot
    self.text += " Modelled %d profiles\n" % ngood
    self.text += " %d profiles could not be created\n" % nbad


  def __str__(self):
    ''' Return as a string. '''
    return self.text

def modeller_statistics(data,
                        experiments,
                        profile_model):
  ''' Return some simple statistics from profile model. '''
  tables = data.split_by_experiment_id()
  assert(len(tables) == len(experiments))
  assert(len(tables) == len(profile_model))
  modeller = profile_model.profiles()
  assert(len(tables) == len(modeller))
  summaries = []
  for index, (table, experiment) in enumerate(zip(tables, experiments)):
    summaries.append(ModellerSummary(
      index,
      table,
      experiment,
      modeller[index]))
  return summaries
