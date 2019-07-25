from __future__ import absolute_import, division, print_function

from collections import namedtuple

from libtbx.phil import parse
from dials.algorithms.refinement import DialsRefineConfigError
from scitbx.array_family import flex
from scitbx import sparse

from dials.algorithms.refinement.restraints.restraints import SingleUnitCellTie
from dials.algorithms.refinement.restraints.restraints import MeanUnitCellTie
from dials.algorithms.refinement.restraints.restraints import LowMemoryMeanUnitCellTie
from dials.algorithms.refinement.restraints.restraints import MedianUnitCellTie

# PHIL options for unit cell restraints
uc_phil_str = """
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
      .help = "Select only the specified experiments when looking up which"
              "parameterisations to apply these restraints to. If an identified"
              "parameterisation affects multiple experiments then the index"
              "of any one of those experiments suffices to restrain that"
              "parameterisation. If None (the default) then the restraints"
              "will be applied to all experiments."
      .type = ints(value_min=0)
  }

  tie_to_group
    .multiple = True
  {
    target = *mean low_memory_mean median
      .type = choice
      .help = "Function to tie group parameter values to"

    sigmas = None
      .help = "The unit cell parameters are associated with sigmas which are"
              "used to determine the weight of each restraint. A sigma of zero"
              "will remove the restraint at that position."
      .type = floats(size=6, value_min=0.)

    id = None
      .help = "Select only the specified experiments when looking up which "
              "parameterisations to apply these restraints to. For every"
              "parameterisation that requires a restraint at least one"
              "experiment index must be supplied. If None (the default) the"
              "restraints will be applied to all experiments."
      .type = ints(value_min=0)
  }
}

"""

uc_phil_scope = parse(uc_phil_str)

# Define a couple of namedtuple types we will use for convenience
ParamIndex = namedtuple("ParamIndex", ["parameterisation", "istart"])
RestraintIndex = namedtuple("RestraintIndex", ["restraint", "istart"])


class RestraintsParameterisation(object):
    def __init__(
        self,
        detector_parameterisations=None,
        beam_parameterisations=None,
        xl_orientation_parameterisations=None,
        xl_unit_cell_parameterisations=None,
        goniometer_parameterisations=None,
    ):

        if detector_parameterisations is None:
            detector_parameterisations = []
        if beam_parameterisations is None:
            beam_parameterisations = []
        if xl_orientation_parameterisations is None:
            xl_orientation_parameterisations = []
        if xl_unit_cell_parameterisations is None:
            xl_unit_cell_parameterisations = []
        if goniometer_parameterisations is None:
            goniometer_parameterisations = []

        # Keep references to all parameterised models
        self._detector_parameterisations = detector_parameterisations
        self._beam_parameterisations = beam_parameterisations
        self._xl_orientation_parameterisations = xl_orientation_parameterisations
        self._xl_unit_cell_parameterisations = xl_unit_cell_parameterisations
        self._goniometer_parameterisations = goniometer_parameterisations

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

        self._exp_to_gon_param = {}
        for gonp in self._goniometer_parameterisations:
            for iexp in gonp.get_experiment_ids():
                self._exp_to_gon_param[iexp] = ParamIndex(gonp, iparam)
            iparam += gonp.num_free()

        # the number of free parameters
        self._nparam = iparam

        # keep a set that will ensure every model parameterisation only gets
        # a single restraint.
        self._param_to_restraint = set()

        # keep lists of restraint objects that we will add
        self._single_model_restraints = []
        self._group_model_restraints = []

    def add_restraints_to_target_xl_unit_cell(self, experiment_id, values, sigma):
        # On input we will have one id value, 6 target values and 6 sigmas.

        # select the right parameterisation, if one exists
        try:
            param_i = self._exp_to_xluc_param[experiment_id]
        except KeyError:
            return

        # fail now if this is already restrained.
        if param_i.parameterisation in self._param_to_restraint:
            raise DialsRefineConfigError(
                "Parameterisation already restrained. Cannot create "
                "additional restraint with experiment {}".format(experiment_id)
            )

        # create new restraint
        tie = SingleUnitCellTie(
            model_parameterisation=param_i.parameterisation, target=values, sigma=sigma
        )

        # add to the restraint list along with the global parameter index
        self._single_model_restraints.append(RestraintIndex(tie, param_i.istart))

        # also add the parameterisation to the set for uniqueness testing
        self._param_to_restraint.add(param_i.parameterisation)

    def add_restraints_to_group_xl_unit_cell(self, target, experiment_ids, sigma):

        # select the right parameterisations, if they exist
        if experiment_ids == "all":
            param_indices = list(self._exp_to_xluc_param.values())
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
                raise DialsRefineConfigError(
                    "Parameterisation already restrained. Cannot create "
                    "additional group restraint for experiment(s) {}".format(
                        str(experiment_ids)
                    )
                )

        # create new group of restraints
        if target == "mean":
            tie = MeanUnitCellTie(model_parameterisations=params, sigma=sigma)
        elif target == "low_memory_mean":
            tie = LowMemoryMeanUnitCellTie(model_parameterisations=params, sigma=sigma)
        elif target == "median":
            tie = MedianUnitCellTie(model_parameterisations=params, sigma=sigma)
        else:
            raise DialsRefineConfigError("target type {} not available".format(target))

        # add to the restraint list along with the global parameter indices
        self._group_model_restraints.append(RestraintIndex(tie, istarts))

    @property
    def num_residuals(self):
        """Get the total number of residuals across all parameterised restraints"""
        n_single = sum(e.restraint.num_residuals for e in self._single_model_restraints)
        n_group = sum(e.restraint.num_residuals for e in self._group_model_restraints)
        return n_single + n_group

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

        # keep track of the row at the start of group models
        group_model_irow = irow

        # process restraints residuals and weights for groups of models
        for r in self._group_model_restraints:
            residuals.extend(flex.double(r.restraint.residuals()))
            weights.extend(flex.double(r.restraint.weights()))

        # now it is clear how many residuals there are we can set up a sparse
        # matrix for the restraints jacobian
        nrows = len(residuals)
        gradients = sparse.matrix(nrows, self._nparam)

        # assign gradients in blocks for the single model restraints
        for irow, r in zip(row_start, self._single_model_restraints):
            icol = r.istart
            # convert square list-of-lists into a 2D array for block assignment
            grads = flex.double(r.restraint.gradients())
            gradients.assign_block(grads, irow, icol)

        # assign gradients in blocks for the group model restraints
        for r in self._group_model_restraints:
            # loop over the included unit cell models, k
            for k, (icol, grads) in enumerate(zip(r.istart, r.restraint.gradients())):
                irow = group_model_irow
                for grad in grads:
                    gradients.assign_block(grad, irow, icol)
                    irow += grad.n_rows
            group_model_irow = irow

        return residuals, gradients, weights
