from __future__ import annotations

import logging

logger = logging.getLogger(__name__)

from libtbx.phil import parse
from scitbx.array_family import flex

from dials.algorithms.refinement import DialsRefineConfigError

phil_str = """
      min_nref_per_parameter = 5
        .help = "the smallest number of reflections per parameter for a"
                "model parameterisation below which the parameterisation will"
                "not be made in full, but the action described below will be"
                "triggered."
        .type = int(value_min=1)

      action = *fail fix remove
        .help = "action to take if there are too few reflections across the"
                "experiments related to a particular model parameterisation."
                "If fail, an exception will be raised and refinement will not"
                "proceed. If fix, refinement will continue but with the"
                "parameters relating to that model remaining fixed at their"
                "initial values. If remove, parameters relating to that model"
                "will be fixed, and in addition all reflections related to"
                "that parameterisation will be removed. This will therefore"
                "remove these reflections from other parameterisations of the"
                "global model too. For example, if a crystal model could not"
                "be parameterised it will be excised completely and not"
                "contribute to the joint refinement of the detector and beam."
                "In the fix mode, reflections emanating from that crystal will"
                "still form residuals and will contribute to detector and beam"
                "refinement."
        .type = choice
"""
phil_scope = parse(phil_str)

# A callback for PredictionParameterisation.get_gradients which will give
# the positions of reflections associated with a particular parameter
def id_associated_refs(result):

    # There are usually 3 parts to results: gradients in X, Y and Z
    vals = list(result.values())
    try:
        vals = [v.as_dense_vector() for v in vals]
    except AttributeError:
        pass

    # Find all non-zero gradients
    vals = [flex.abs(v) > 1e-10 for v in vals]

    # Find any reflection that contributed a non-zero gradient
    val = vals[0]
    for v in vals[1:]:
        val = val | v
    return val


# A callback for PredictionParameterisation.get_gradients which will count
# the number of reflections associated with a particular parameter
def count_associated_refs(result):
    refs = id_associated_refs(result)
    return refs.count(True)


class AutoReduce:
    """Checks for over-parameterisation of models and acts in that case.

    Tests each model parameterisation to ensure there are enough
    reflections in refinement to support that parameterisation. If there are
    not then some action is taken. More details are given in documentation
    within the phil_str alongside this class definition.
    """

    def __init__(
        self,
        options,
        pred_param,
        reflection_manager,
        constraints_manager=None,
        constraints_manager_factory=None,
    ):
        """Initialise the AutoReduce object

        Args:
            options: A PHIL scope containing the auto reduction options
            pred_param: A PredictionParameterisation object containing the
                global model parameterisation
            reflection_manager: A ReflectionManager object tracking the observations
                for use in refinement
            constraints_manager: An optional ConstraintsManager object to allow correct
                determination of the number of reflections per parameter when constraints
                are present
            constraints_manager_factory: An optional ConstraintsManagerFactory, required
                if constraints_manager is present, to recreate the constraints_manager
                if parameters have been fixed
        """

        self._options = options
        self.pred_param = pred_param
        self.param_names = self.pred_param.get_param_names()
        self.reflection_manager = reflection_manager
        self.constraints_manager = constraints_manager
        self.constraints_manager_factory = constraints_manager_factory

    def _nref_per_param(self):
        # Calculate the number of reflections for each parameter, taking
        # constraints into account
        obs = self.reflection_manager.get_obs()
        try:
            self.pred_param.compose(obs)
        except AttributeError:
            pass

        nref_per_param = flex.size_t(
            self.pred_param.get_gradients(obs, callback=count_associated_refs)
        )

        if self.constraints_manager is not None:
            for link in self.constraints_manager.get_constrained_parameter_indices():
                sel = flex.size_t(link)
                total = flex.sum(nref_per_param.select(sel))
                nref_per_param.set_selected(sel, total)

        return nref_per_param

    def check_and_fail(self):
        """Check for too few reflections to support the model parameterisation.

        Test each parameterisation of each type against the reflections it affects.

        Returns:
            None

        Raises:
            DialsRefineConfigError: If there are too few reflections to support
            a parameterisation.
        """

        sel = (
            self._nref_per_param() < self._options.min_nref_per_parameter
        ).iselection()
        if len(sel) > 0:
            names = ", ".join([self.param_names[i] for i in sel])
            msg = f"Too few reflections to parameterise {names}.\n"
            msg += (
                "Try setting "
                "refinement.parameterisation.auto_reduction.action "
                "to fix these parameters (=fix) or additionally remove the "
                "associated reflections (=remove)."
            )
            raise DialsRefineConfigError(msg)

    def check_and_fix(self):
        """Fix parameters when there are too few reflections.

        Test each parameter against the reflections it affects and fix any for
        which there are too few reflections.

        Returns:
            None
        """

        sel = self._nref_per_param() < self._options.min_nref_per_parameter
        isel = sel.iselection()
        if len(isel) > 0:
            names = ", ".join([self.param_names[i] for i in isel])
            msg = f"Too few reflections to parameterise {names}.\n"
            msg += "These parameters will be fixed for refinement."
            logger.warning(msg)
        self.pred_param.fix_params(sel)
        if self.constraints_manager is not None:
            self.constraints_manager = self.constraints_manager_factory()

    def check_and_remove(self):
        """Fix parameters and remove reflections when there are too few reflections.

        Test each parameter against the reflections it affects and fix any for
        which there are too few reflections. In addition, remove all reflections
        that are associated with that parameter to ensure they play no part in
        refinement. This process is iterative.

        Returns:
            None

        Raises:
            DialsRefineConfigError: error if only one single panel detector is present.
        """

        # If there is only one detector in a single experiment, the detector should
        # be multi-panel for remove to make sense
        det_params = self.pred_param.get_detector_parameterisations()
        if len(det_params) == 1:
            n_exp = len(det_params[0].get_experiment_ids())
            if n_exp == 1 and not det_params[0].is_multi_state():
                raise DialsRefineConfigError(
                    "For single experiment, single panel refinement "
                    "auto_reduction.action=remove cannot be used as it could only "
                    "remove all reflections from refinement"
                )

        while True:
            obs = self.reflection_manager.get_obs()
            try:
                self.pred_param.compose(obs)
            except AttributeError:
                pass

            refs_by_parameters = self.pred_param.get_gradients(
                obs, callback=id_associated_refs
            )
            nref_per_param = flex.size_t(
                [refs.count(True) for refs in refs_by_parameters]
            )

            if self.constraints_manager is not None:
                for (
                    link
                ) in self.constraints_manager.get_constrained_parameter_indices():
                    sel = flex.size_t(link)
                    total = flex.sum(nref_per_param.select(sel))
                    nref_per_param.set_selected(sel, total)

            sel = nref_per_param < self._options.min_nref_per_parameter
            if sel.count(True) == 0:
                break

            names = ", ".join([self.param_names[i] for i in sel.iselection()])
            msg = f"Too few reflections to parameterise {names}.\n"
            msg += (
                "These parameters will be fixed for refinement and "
                "the associated reflections will be removed."
            )
            logger.warning(msg)

            self.pred_param.fix_params(sel)

            if self.constraints_manager is not None:
                self.constraints_manager = self.constraints_manager_factory()

            refs_to_filter = flex.bool(len(obs), True)
            for remove, refs in zip(sel, refs_by_parameters):
                if remove:
                    refs_to_filter = refs_to_filter & ~refs

            # only keep refs not associated with this parameterisation
            self.reflection_manager.filter_obs(refs_to_filter)

    def __call__(self):
        """Perform checks and parameter reduction according to the selected option.

        Returns:
            None
        """

        if self._options.action == "fail":
            self.check_and_fail()

        elif self._options.action == "fix":
            self.check_and_fix()

        elif self._options.action == "remove":
            self.check_and_remove()
