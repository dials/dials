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
from libtbx.utils import Sorry
from scitbx.array_family import flex
from scitbx import sparse

from dials.algorithms.refinement.restraints.restraints import SingleUnitCellTie

# PHIL options for unit cell restraints
uc_phil_str = '''
restraints
  .help = "Least squares restraints to use in refinement."
  .expert_level = 1
{
  tie_to_target
    .multiple = True
  {
    values = None
      .type = floats(size=6)
      .help = "Target unit cell parameters for the restraint for this"
              "parameterisation"

    sigmas = None
      .help = "The unit cell target values are associated with sigmas which are"
              "used to determine the weight of each restraint. A sigma of zero"
              "or None will remove the restraint at that position"
      .type = floats(size=6, value_min=0.)

    id = None
      .help = "Index of an experiment affected by this restraint. If the"
              "relevant parameterisation affects multiple experiments then any"
              "one of the indices may be supplied"
      .type = int(value_min=0)
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
      .help = "Indices of experiments affected by this restraint. For every"
              "parameterisation that requires a restraint at least one"
              "experiment index must be supplied."
      .type = ints(value_min=0)
  }
}

'''

uc_phil_scope = parse(uc_phil_str)

# Define a couple of namedtuple types we will use for convenience
from collections import namedtuple
ParamIndex = namedtuple('ParamIndex', ['parameterisation', 'istart'])
RestraintIndex = namedtuple('RestraintIndex', ['restraint', 'istart'])

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
    # global parameter index for each that tells us which parameters have
    # non-zero derivatives
    iparam = 0

    self._exp_to_det_param = {}
    for detp in self._detector_parameterisations:
      for iexp in detp.get_experiment_ids():
        self._exp_to_det_param[iexp] = ParamIndex(detp, iparam)
      iparam += detp.num_free()

    self._exp_to_beam_param = {}
    for beamp in self._beam_parameterisations:
      for iexp in beamp.get_experiment_ids():
        self._exp_to_beam_param[iexp] = ParamIndex(beamp, iparam)
      iparam += beamp.num_free()

    self._exp_to_xlo_param = {}
    for xlop in self._xl_orientation_parameterisations:
      for iexp in xlop.get_experiment_ids():
        self._exp_to_xlo_param[iexp] = ParamIndex(xlop, iparam)
      iparam += xlop.num_free()

    self._exp_to_xluc_param = {}
    for xlucp in self._xl_unit_cell_parameterisations:
      for iexp in xlucp.get_experiment_ids():
        self._exp_to_xluc_param[iexp] = ParamIndex(xlucp, iparam)
      iparam += xlucp.num_free()

    # the number of free parameters
    self._nparam = iparam

    # keep a set that will ensure every model parameterisation only gets
    # a single restraint.
    self._param_to_restraint = set()

    # keep a list of restraint objects that we will add
    self._restraints = []

    return

  #def add_restraints_to_target_detector(self):
  #
  #  return
  #
  #def add_restraints_to_target_beam(self):
  #
  #  return
  #
  #def add_restraints_to_target_xl_orientation(self):
  #
  #  return

  def add_restraints_to_target_xl_unit_cell(self,
    experiment_id, values, sigma):

    # On input we will have one id value, 6 target values and 6 sigmas.

    # select the right parameterisation
    param_i = self._exp_to_xluc_param[experiment_id]

    # fail now if this is already restrained.
    if param_i.parameterisation in self._param_to_restraint:
      raise Sorry("Parameterisation already restrained. Cannot create "
                  "additional restraint with experiment {0}".format(experiment_id))

    # create new restraint
    tie = SingleUnitCellTie(model_parameterisation=param_i.parameterisation,
                            target=values,
                            sigma=sigma)

    # add to the restraint list along with the global parameter index
    self._restraints.append(RestraintIndex(tie, param_i.istart))

    # also add the parameterisation to the set for uniqueness testing
    self._param_to_restraint.add(param_i.parameterisation)

    return

  #def add_restraints_to_group_detector(self):
  #
  #  return
  #
  #def add_restraints_to_group_beam(self):
  #
  #  return
  #
  #def add_restraints_to_group_xl_orientation(self):
  #
  #  return

  def add_restraints_to_group_xl_unit_cell(self):

    pass
    return

  def get_values_and_gradients(self):

    values = []
    for r in self._restraints:
      values.extend(r.restraint.values())

    nrows = len(values)
    gradients = sparse.matrix(nrows, self._nparam)

    for irow, r in enumerate(self._restraints):
      icol = r.istart
      # convert square list-of-lists into a 2D array for block assignment
      grads = flex.double(r.restraint.gradients())
      gradients.assign_block(grads, irow, icol)

    return values, gradients
