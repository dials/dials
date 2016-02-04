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
from dials.algorithms.refinement.restraints.restraints import MeanUnitCellTie

# PHIL options for unit cell restraints
uc_phil_str = '''
restraints
  .help = "Least squares unit cell restraints to use in refinement."
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
              "will remove the restraint at that position. If symmetry"
              "constrains two cell dimensions to be equal then only the"
              "smaller of the two sigmas will be kept"
      .type = floats(size=6, value_min=0.)

    id = None
      .help = "Index of experiments affected by this restraint to look up which"
              "parameterisations to apply the restraint to. If an identified"
              "parameterisation affects multiple experiments then the index"
              "of any one of those experiments suffices to restrain that"
              "parameterisation."
      .type = ints(value_min=0)
  }

  tie_to_group
    .multiple = True
  {
    target = *mean median
      .type = choice
      .help = "Function to tie group parameter values to"

    sigmas = None
      .help = "The unit cell parameters are associated with sigmas which are"
              "used to determine the weight of each restraint. A sigma of zero"
              "will remove the restraint at that position."
      .type = floats(size=6, value_min=0.)

    id = None
      .help = "Indices of experiments affected by this restraint. For every"
              "parameterisation that requires a restraint at least one"
              "experiment index must be supplied, unless using apply_to_all"
      .type = ints(value_min=0)

    apply_to_all = False
      .help = "Shorthand to restrain the unit cells across all experiments"
      .type = bool
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

    # keep lists of restraint objects that we will add
    self._single_model_restraints = []
    self._group_model_restraints = []

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

    # select the right parameterisation, if one exists
    try:
      param_i = self._exp_to_xluc_param[experiment_id]
    except KeyError:
      return

    # fail now if this is already restrained.
    if param_i.parameterisation in self._param_to_restraint:
      raise Sorry("Parameterisation already restrained. Cannot create "
                  "additional restraint with experiment {0}".format(experiment_id))

    # create new restraint
    tie = SingleUnitCellTie(model_parameterisation=param_i.parameterisation,
                            target=values,
                            sigma=sigma)

    # add to the restraint list along with the global parameter index
    self._single_model_restraints.append(RestraintIndex(tie, param_i.istart))

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

  def add_restraints_to_group_xl_unit_cell(self, target, experiment_ids, sigma):

    # select the right parameterisations, if they exist
    if experiment_ids == 'all':
      param_indices = self._exp_to_xluc_param.values()
    else:
      param_indices = []
      for exp_id in experiment_ids:
        try:
          param_indices.append(self._exp_to_xluc_param[exp_id])
        except KeyError:
          # ignore experiment without a parameterisation
          pass
    params = [e.parameterisation for e in param_indices]
    istarts = [e.istart for e in param_indices]

    # fail if any of the parameterisations has already been restrained.
    for param in params:
      if param in self._param_to_restraint:
        raise Sorry("Parameterisation already restrained. Cannot create "
                    "additional group restraint for experiment(s) {0}".format(str(
                      param_i.parameterisation.get_experiment_ids())))

    # create new group of restraints
    if target == 'mean':
      tie = MeanUnitCellTie(model_parameterisations=params,
                            sigma=sigma)
    else:
      raise Sorry("target type {0} not available".format(target))

    # add to the restraint list along with the global parameter indices
    self._group_model_restraints.append(RestraintIndex(tie, istarts))

    return

  def get_residuals_gradients_and_weights(self):

    residuals = flex.double()
    weights = flex.double()
    row_start = []
    irow = 0

    # process restraints residuals and weights for single models
    for r in self._single_model_restraints:
      res = r.restraint.residuals()
      wgt = r.restraint.weights()
      residuals.extend(flex.double(res))
      weights.extend(flex.double(wgt))
      row_start.append(irow)
      irow += len(res)

    group_model_irow = irow

    # process restraints residuals and weights for groups of models
    for r in self._group_model_restraints:
      residuals.extend(flex.double(r.restraint.residuals()))
      weights.extend(flex.double(r.restraint.weights()))

    # set up a sparse matrix for the restraints jacobian
    nrows = len(residuals)
    gradients = sparse.matrix(nrows, self._nparam)

    # assign gradients in blocks for the single model restraints
    for irow, r in zip(row_start, self._single_model_restraints):
      icol = r.istart
      # convert square list-of-lists into a 2D array for block assignment
      grads = flex.double(r.restraint.gradients())
      gradients.assign_block(grads, irow, icol)

    # assign gradients in blocks for the group model restraints
    irow = group_model_irow
    for r in self._group_model_restraints:
      # loop over the included unit cell models, k
      for k, (icol, grads) in enumerate(zip(r.istart, r.restraint.gradients())):
        this_cell_grads = flex.double(grads)
        other_cells_grads = flex.double(r.restraint.gradients_of_the_mean(k)) * -1.0
        # the cell model k has a parameterisation that starts with column icol.
        # the block of gradients of its residuals wrt its parameters is the
        # block this_cell_grads, whilst the block of gradients of other cell's
        # residuals wrt its parameters is other_cells_grads.
        #
        # need to loop over the starting rows and write the blocks
        for j in range(r.restraint._nxls):
          jrow = irow + j * r.restraint.nrestraints_per_cell
          if j == k:
            gradients.assign_block(this_cell_grads, jrow, icol)
          else:
            gradients.assign_block(other_cells_grads, jrow, icol)
      irow += r.restraint._nxls * r.restraint.nrestraints_per_cell

    return residuals, gradients, weights
