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
from libtbx.phil import parse

# The phil scope for restraints
phil_str = '''
restraints
  .help = "Least squares restraints to use in refinement."
  .expert_level = 1
{
  tie_to_target
    .multiple = True
  {
    values = None
      .type = floats()
      .help = "Target values for the restraint for this parameterisation"

    sigmas = None
      .help = "A sigma of zero or None will remove the restraint at that position"
      .type = floats(value_min=0.)

    id = None
      .help = "Experiment indices of experiments to be affected by this restraint"
      .type = ints(value_min=0)
  }

  tie_to_group
    .multiple = True
  {
    target = *mean median
      .type = choice
      .help = "Function to tie group parameter values to"

    sigmas = None
      .help = "A sigma of zero or None will remove the restraint at that position"
      .type = floats(value_min=0.)

    id = None
      .help = "Experiment indices of experiments to be affected by this restraint"
      .type = ints(value_min=0)
  }
}

'''

phil_scope = parse(phil_str)

class RestraintsParameterisation(object):

  def __init__(self, detector_parameterisations = None,
               beam_parameterisations = None,
               xl_orientation_parameterisations = None,
               xl_unit_cell_parameterisations = None):

    # Keep references to all parameterised models
    self._detector_parameterisations = detector_parameterisations
    self._beam_parameterisations = beam_parameterisations
    self._xl_orientation_parameterisations = xl_orientation_parameterisations
    self._xl_unit_cell_parameterisations = xl_unit_cell_parameterisations

    # Loop over all parameterisations, extract experiment IDs and record
    # global parameter index for each
    iparam = 0
    from collections import namedtuple
    ParamIndex = namedtuple('ParamIndex', ['parameterisation', 'istart'])

    exp_to_det_param = {}
    for detp in self._detector_parameterisations:
      for iexp in detp.get_experiment_ids():
        exp_to_det_param[iexp] = ParamIndex(detp, iparam)
      iparam += detp.num_free()

    exp_to_beam_param = {}
    for beamp in self._beam_parameterisations:
      for iexp in beamp.get_experiment_ids():
        exp_to_beam_param[iexp] = ParamIndex(beamp, iparam)
      iparam += beamp.num_free()

    exp_to_xlo_param = {}
    for xlop in self._xl_orientation_parameterisations:
      for iexp in xlop.get_experiment_ids():
        exp_to_xlo_param[iexp] = ParamIndex(xlop, iparam)
      iparam += xlop.num_free()

    exp_to_xluc_param = {}
    for xlucp in self._xl_unit_cell_parameterisations:
      for iexp in xlucp.get_experiment_ids():
        exp_to_xluc_param[iexp] = ParamIndex(xlucp, iparam)
      iparam += xlucp.num_free()

    #from dials.util.command_line import interactive_console; interactive_console()
