#!/usr/bin/env python
#
#  restraints_parameterisation.py
#
#  Copyright (C) 2015 Diamond Light Source and STFC Rutherford Appleton
#                     Laboratory, UK.
#
#  Author: David Waterman
#
#  This code is distributed under the BSD license, a copy of which is
#  included in the root directory of this package.

from __future__ import division

# The phil scope for restraints
phil_str = '''
restraints
  .help = "Least squares restraints to use in refinement."
  .expert_level = 1
{
  tie_to_target
    .help = ""
  {
    values = None
      .type = floats()
      .help = "Target values for the restraint for this parameterisation"

    sigmas = None
      .help = "A sigma of zero or None will remove the restraint at that position"
      .type = floats(value_min=0.)

    experiments = None
      .help = "Experiment indices of experiments to be affected by this restraint"
      .type = ints(value_min=0)
      .
  }

  tie_to_group
    .help = ""
  {
    target = *mean median
      .type = choice
      .help = "Function to tie group parameter values to"

    sigmas = None
      .help = "A sigma of zero or None will remove the restraint at that position"
      .type = floats(value_min=0.)

    experiments = None
      .help = "Experiment indices of experiments to be affected by this restraint"
      .type = ints(value_min=0)
      .
  }
}

'''

class RestraintsParameterisation(object):

  def __init__(self):
