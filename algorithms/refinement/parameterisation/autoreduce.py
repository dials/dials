from __future__ import absolute_import, division, print_function
import logging

logger = logging.getLogger(__name__)

from scitbx.array_family import flex
from dials.algorithms.refinement import DialsRefineConfigError

# PHIL
from libtbx.phil import parse

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

      detector_reduce = False
        .type = bool
        .help = "Special case designed for detector metrology refinement"
                "(particularly of the CSPAD). See detector_reduce_list for"
                "details."
        .expert_level = 2

      detector_reduce_list = Dist Tau2 Tau3
        .type = strings
        .help = "Partial names to match to detector parameters to try fixing."
                "If there are still not"
                "enough parameters for refinement after fixing these, then"
                "fail. This is to ensure that metrology refinement never"
                "completes if it is not able to refine some panels. The default"
                "is to try fixing the distance as well as Tau2 and Tau3"
                "rotations of detector panel, leaving the in-plane shifts and"
                "the rotation around the detector normal for refinement."
                "groups only."
        .expert_level = 2
"""
phil_scope = parse(phil_str)

# A callback for PredictionParameterisation.get_gradients which will give
# the positions of reflections associated with a particular parameter
def id_associated_refs(result):

    # There are usually 3 parts to results: gradients in X, Y and Z
    vals = [v for v in result.values()]
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


class AutoReduce(object):
    """Checks for over-parameterisation of models and acts in that case.

    Tests each model parameterisation to ensure there are enough
    reflections in refinement to support that parameterisation. If there are
    not then some action is taken. More details are given in documentation
    within the phil_str alongside this class definition.

    Attributes:
        det_params (list): A list of DetectorParameterisation objects
        beam_params (list): A list of BeamParameterisation objects
        xl_ori_params (list): A list of CrystalOrientationParameterisation objects
        xl_uc_params (list): A list of CrystalUnitCellParameterisation objects
        gon_params (list): A list of GoniometerParameterisation objects
        reflections: A reflection table
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
            det_params (list): A list of DetectorParameterisation objects
            beam_params (list): A list of BeamParameterisation objects
            xl_ori_params (list): A list of CrystalOrientationParameterisation
                objects
            xl_uc_params (list): A list of CrystalUnitCellParameterisation objects
            gon_params (list): A list of GoniometerParameterisation objects
            reflection_manager: The ReflectionManager object handling reflection
                data for refinement
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
            msg = "Too few reflections to parameterise {0}.\n".format(names)
            msg += "Try modifying refinement.parameterisation.auto_reduction options"
            raise DialsRefineConfigError(msg)

    def check_and_fix(self):
        """Fix parameters when there are too few reflections.

        Test each parameterisation of each type against the reflections it affects.
        If there are too few reflections to support that parameterisation, fix the
        parameters.

        Returns:
            None
        """

        sel = self._nref_per_param() < self._options.min_nref_per_parameter
        isel = sel.iselection()
        if len(isel) > 0:
            names = ", ".join([self.param_names[i] for i in isel])
            msg = "Too few reflections to parameterise {0}.\n".format(names)
            msg += "These parameters will be fixed for refinement."
            logger.warning(msg)
        self.pred_param.fix_params(sel)
        if self.constraints_manager is not None:
            self.constraints_manager = self.constraints_manager_factory()

    def check_and_remove(self):
        """Fix parameters and remove reflections when there are too few reflections.

        Test each parameterisation of each type against the reflections it affects.
        If there are too few reflections to support that parameterisation, fix the
        parameters and remove those reflections so that they will not be included
        in refinement.

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

            refs_by_parameterisation = self.pred_param.get_gradients(
                obs, callback=id_associated_refs
            )
            nref_per_param = flex.size_t(
                [refs.count(True) for refs in refs_by_parameterisation]
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
            msg = "Too few reflections to parameterise {0}.\n".format(names)
            msg += (
                "These parameters will be fixed for refinement and "
                "the associated reflections will be removed."
            )
            logger.warning(msg)

            self.pred_param.fix_params(sel)

            if self.constraints_manager is not None:
                self.constraints_manager = self.constraints_manager_factory()

            for remove, refs in zip(sel, refs_by_parameterisation):
                if remove:
                    # only keep refs not associated with this parameterisation
                    self.reflection_manager.filter_obs(~refs)

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
