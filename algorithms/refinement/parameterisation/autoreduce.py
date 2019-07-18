from __future__ import absolute_import, division, print_function
import logging

logger = logging.getLogger(__name__)

# Parameterisation auto reduction helpers
from dials_refinement_helpers_ext import surpl_iter as surpl
from dials_refinement_helpers_ext import uc_surpl_iter as uc_surpl
from dials_refinement_helpers_ext import pg_surpl_iter as pg_surpl

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


class AutoReduce(object):
    """Checks for over-parameterisation of models and acts in that case.

    Tests each provided model parameterisation to ensure there are enough
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
        det_params,
        beam_params,
        xl_ori_params,
        xl_uc_params,
        gon_params,
        reflection_manager,
        scan_varying=False,
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
            scan_varying (bool): Whether preparing for scan-varying refinement or
                scan static refinement

        """

        self.det_params = det_params
        self.beam_params = beam_params
        self.xl_ori_params = xl_ori_params
        self.xl_uc_params = xl_uc_params
        self.gon_params = gon_params
        self.reflection_manager = reflection_manager

        self._options = options
        self._scan_varying = scan_varying

        # A template logging message to fill in when failing
        self._failmsg = (
            "Too few reflections to parameterise {0}\nTry modifying "
            "refinement.parameterisation.auto_reduction options"
        )

    # Determine if there are enough reflections to support a particular
    # parameterisation. First, a minimum number of reflections is determined,
    # by the product of the number of free parameters and a user-provided
    # minimum number of reflections per parameter. The total number of reflections
    # affected by this parameterisation is calculated, and the difference between
    # that and the minimum number of reflections is returned.
    def _surplus_reflections(self, p):
        reflections = self.reflection_manager.get_obs()
        cutoff = self._options.min_nref_per_parameter * p.num_free()
        return surpl(reflections["id"], p.get_experiment_ids()).result - cutoff

    # Special version of _surplus_reflections for crystal unit cell
    # parameterisations. In some cases certain parameters of a unit cell
    # parameterisation may affect only some subset of the total number of
    # reflections. For example, for an orthorhombic cell the g_param_0 parameter
    # has no effect on predictions in the plane (0,k,l). Here, take the number
    # of affected reflections for each parameter individually into account.
    def _unit_cell_surplus_reflections(self, p):
        F_dbdp = flex.mat3_double(p.get_ds_dp())
        min_nref = self._options.min_nref_per_parameter
        reflections = self.reflection_manager.get_obs()
        # if no free parameters, do as _surplus_reflections
        if len(F_dbdp) == 0:
            exp_ids = p.get_experiment_ids()
            isel = flex.size_t()
            for exp_id in exp_ids:
                isel.extend((reflections["id"] == exp_id).iselection())
            return len(isel)
        return (
            uc_surpl(
                reflections["id"],
                reflections["miller_index"],
                p.get_experiment_ids(),
                F_dbdp,
            ).result
            - min_nref
        )

    # Special version of _surplus_reflections for hierarchical multi-panel detector
    # parameterisations. In that case, certain parameters affect only the
    # reflections that fall on a particular panel group of the detector.
    def _panel_gp_surplus_reflections(self, p, pnl_ids, group):
        exp_ids = p.get_experiment_ids()
        gp_params = [gp == group for gp in p.get_param_panel_groups()]
        fixlist = p.get_fixed()
        free_gp_params = [a and not b for a, b in zip(gp_params, fixlist)]
        nparam = free_gp_params.count(True)
        cutoff = self._options.min_nref_per_parameter * nparam
        reflections = self.reflection_manager.get_obs()
        surplus = pg_surpl(
            reflections["id"], reflections["panel"], pnl_ids, exp_ids, cutoff
        ).result
        return surplus

    def _weak_parameterisation_search(self):
        weak = None
        nref_deficit = 0
        panels = None
        pnl_gp = None
        name = None
        for i, p in enumerate(self.beam_params):
            net_nref = self._surplus_reflections(p)
            if net_nref < nref_deficit:
                nref_deficit = net_nref
                weak = p
                name = "Beam{}".format(i + 1)
        for i, p in enumerate(self.xl_ori_params):
            net_nref = self._surplus_reflections(p)
            if net_nref < nref_deficit:
                nref_deficit = net_nref
                weak = p
                name = "Crystal{} orientation".format(i + 1)
        for i, p in enumerate(self.xl_uc_params):
            net_nref = self._unit_cell_surplus_reflections(p)
            if net_nref < nref_deficit:
                nref_deficit = net_nref
                weak = p
                name = "Crystal{} unit cell".format(i + 1)
        for i, p in enumerate(self.det_params):
            try:
                pnl_groups = p.get_panel_ids_by_group()
                for igp, gp in enumerate(pnl_groups):
                    net_nref = self._panel_gp_surplus_reflections(p, gp, igp)
                    if net_nref < nref_deficit:
                        nref_deficit = net_nref
                        weak = p
                        panels = gp
                        pnl_gp = igp
                        name = "Detector{0}PanelGroup{1}".format(i + 1, pnl_gp + 1)
            except AttributeError:  # non-hierarchical detector parameterisation
                net_nref = self._surplus_reflections(p)
                if net_nref < nref_deficit:
                    nref_deficit = net_nref
                    weak = p
                    panels = None
                    pnl_gp = None
                    name = "Detector{}".format(i + 1)
        for i, p in enumerate(self.gon_params):
            net_nref = self._surplus_reflections(p)
            if net_nref < nref_deficit:
                nref_deficit = net_nref
                weak = p
                name = "Goniometer{}".format(i + 1)
        return {
            "parameterisation": weak,
            "panels": panels,
            "panel_group_id": pnl_gp,
            "name": name,
        }

    def detector_reduce(self):
        """Reduce detector parameters.

        Special case intended for metrology refinement of multi-panel detectors."""
        reduce_list = self._options.detector_reduce_list
        for i, dp in enumerate(self.det_params):
            to_fix = flex.bool(dp.get_fixed())
            try:  # test for hierarchical detector parameterisation
                pnl_groups = dp.get_panel_ids_by_group()
                for igp, gp in enumerate(pnl_groups):
                    surplus = self._panel_gp_surplus_reflections(dp, gp, igp)
                    if surplus < 0:
                        msg = (
                            "Require {0} more reflections to parameterise Detector{1} "
                            "panel group {2}"
                        )
                        logger.warning(
                            msg.format(-1 * surplus, i + 1, igp + 1)
                            + "\nAttempting reduction of non-essential parameters"
                        )
                        names = cls._filter_parameter_names(dp)
                        prefix = "Group{}".format(igp + 1)
                        reduce_this_group = [prefix + e for e in reduce_list]
                        to_fix |= flex.bool(string_sel(reduce_this_group, names))
                        # try again, and fail if still unsuccessful
                        surplus = self._panel_gp_surplus_reflections(dp, gp, igp)
                        if surplus < 0:
                            msg = msg.format(-1 * surplus, i + 1, igp + 1)
                            raise DialsRefineConfigError(msg + "\nFailing.")
            except AttributeError:
                if self._surplus_reflections(dp) < 0:
                    mdl = "Detector{}".format(i + 1)
                    msg = self._failmsg.format(mdl)
                    raise DialsRefineConfigError(msg)
            dp.set_fixed(to_fix)

    def check_and_fail(self):
        """Check for too few reflections to support the model parameterisation.

        Test each parameterisation of each type against the reflections it affects.

        Returns:
            None

        Raises:
            DialsRefineConfigError: If there are too few reflections to support
            a parameterisation.
        """

        for i, bp in enumerate(self.beam_params):
            if self._surplus_reflections(bp) < 0:
                mdl = "Beam{}".format(i + 1)
                msg = self._failmsg.format(mdl)
                raise DialsRefineConfigError(msg)

        for i, xlo in enumerate(self.xl_ori_params):
            if self._surplus_reflections(xlo) < 0:
                mdl = "Crystal{} orientation".format(i + 1)
                msg = self._failmsg.format(mdl)
                raise DialsRefineConfigError(msg)

        for i, xluc in enumerate(self.xl_uc_params):
            if self._unit_cell_surplus_reflections(xluc) < 0:
                mdl = "Crystal{} unit cell".format(i + 1)
                msg = self._failmsg.format(mdl)
                raise DialsRefineConfigError(msg)

        for i, dp in enumerate(self.det_params):
            try:  # test for hierarchical detector parameterisation
                pnl_groups = dp.get_panel_ids_by_group()
                for igp, gp in enumerate(pnl_groups):
                    if self._panel_gp_surplus_reflections(dp, gp, igp) < 0:
                        msg = "Too few reflections to parameterise Detector{0} panel group {1}"
                        msg = msg.format(i + 1, igp + 1)
                        msg += "\nTry modifying refinement.parameterisation.auto_reduction options"
                        raise DialsRefineConfigError(msg)
            except AttributeError:
                if self._surplus_reflections(dp) < 0:
                    mdl = "Detector{}".format(i + 1)
                    msg = self._failmsg.format(mdl)
                    raise DialsRefineConfigError(msg)

        for i, gonp in enumerate(self.gon_params):
            if self._surplus_reflections(gonp) < 0:
                mdl = "Goniometer{}".format(i + 1)
                msg = self._failmsg.format(mdl)
                raise DialsRefineConfigError(msg)

    def check_and_fix(self):
        """Fix parameters when there are too few reflections.

        Test each parameterisation of each type against the reflections it affects.
        If there are too few reflections to support that parameterisation, fix the
        parameters.

        Returns:
            None
        """
        warnmsg = "Too few reflections to parameterise {0}"
        tmp = []
        for i, bp in enumerate(self.beam_params):
            if self._surplus_reflections(bp) >= 0:
                tmp.append(bp)
            else:
                mdl = "Beam{}".format(i + 1)
                msg = warnmsg.format(mdl)
                logger.warning(msg)
        self.beam_params = tmp

        tmp = []
        for i, xlo in enumerate(self.xl_ori_params):
            if self._surplus_reflections(xlo) >= 0:
                tmp.append(xlo)
            else:
                mdl = "Crystal{} orientation".format(i + 1)
                msg = warnmsg.format(mdl)
                logger.warning(msg)
        self.xl_ori_params = tmp

        tmp = []
        for i, xluc in enumerate(self.xl_uc_params):
            if self._unit_cell_surplus_reflections(xluc) >= 0:
                tmp.append(xluc)
            else:
                mdl = "Crystal{} unit cell".format(i + 1)
                msg = warnmsg.format(mdl)
                logger.warning(msg)
        self.xl_uc_params = tmp

        tmp = []
        for i, dp in enumerate(self.det_params):
            fixlist = dp.get_fixed()
            try:  # test for hierarchical detector parameterisation
                pnl_groups = dp.get_panel_ids_by_group()
                for igp, gp in enumerate(pnl_groups):
                    if self._panel_gp_surplus_reflections(dp, gp, igp) < 0:
                        msg = "Too few reflections to parameterise Detector{0}PanelGroup{1}"
                        msg = msg.format(i + 1, igp + 1)
                        logger.warning(msg)
                        gp_params = [gp == igp for gp in dp.get_param_panel_groups()]
                        for j, val in enumerate(gp_params):
                            if val:
                                fixlist[j] = True
                dp.set_fixed(fixlist)
                if dp.num_free() > 0:
                    tmp.append(dp)
                else:
                    msg = "No parameters remain free for Detector{}".format(i + 1)
                    logger.warning(msg)
            except AttributeError:
                if self._surplus_reflections(dp) >= 0:
                    tmp.append(dp)
                else:
                    mdl = "Detector{}".format(i + 1)
                    msg = warnmsg.format(mdl)
                    logger.warning(msg)
        self.det_params = tmp

        tmp = []
        for i, gonp in enumerate(self.gon_params):
            if self._surplus_reflections(gonp) >= 0:
                tmp.append(gonp)
            else:
                mdl = "Goniometer{}".format(i + 1)
                msg = warnmsg.format(mdl)
                logger.warning(msg)
        self.gon_params = tmp

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
        if len(self.det_params) == 1:
            n_exp = len(self.det_params[0].get_experiment_ids())
            if n_exp == 1 and not self.det_params[0].is_multi_state():
                raise DialsRefineConfigError(
                    "For single experiment, single panel refinement "
                    "auto_reduction.action=remove cannot be used as it could only "
                    "remove all reflections from refinement"
                )

        # Define a warning message template to use each search iteration
        warnmsg = "Too few reflections to parameterise {0}"
        warnmsg += (
            "\nAssociated reflections will be removed from the Reflection Manager"
        )

        while True:
            # Identify a poorly-supported parameterisation
            dat = self._weak_parameterisation_search()
            if dat["parameterisation"] is None:
                break
            exp_ids = dat["parameterisation"].get_experiment_ids()
            msg = warnmsg.format(dat["name"])

            # Fix relevant parameters and identify observations to remove
            obs = self.reflection_manager.get_obs()
            isel = flex.size_t()
            if dat["panels"] is not None:
                fixlist = dat["parameterisation"].get_fixed()
                pnl_gps = dat["parameterisation"].get_param_panel_groups()
                for i, gp in enumerate(pnl_gps):
                    if gp == dat["panel_group_id"]:
                        fixlist[i] = True
                dat["parameterisation"].set_fixed(fixlist)
                # identify observations on this panel group from associated experiments
                for exp_id in exp_ids:
                    subsel = (obs["id"] == exp_id).iselection()
                    panels_this_exp = obs["panel"].select(subsel)
                    for pnl in dat["panels"]:
                        isel.extend(subsel.select(panels_this_exp == pnl))
            else:
                fixlist = [True] * dat["parameterisation"].num_total()
                dat["parameterisation"].set_fixed(fixlist)
                # identify observations from the associated experiments
                for exp_id in exp_ids:
                    isel.extend((obs["id"] == exp_id).iselection())

            # Now remove the selected reflections
            sel = flex.bool(len(obs), True)
            sel.set_selected(isel, False)
            self.reflection_manager.filter_obs(sel)
            logger.warning(msg)

        # Strip out parameterisations with zero free parameters
        self.beam_params = [p for p in self.beam_params if p.num_free() > 0]
        self.xl_ori_params = [p for p in self.xl_ori_params if p.num_free() > 0]
        self.xl_uc_params = [p for p in self.xl_uc_params if p.num_free() > 0]
        self.det_params = [p for p in self.det_params if p.num_free() > 0]
        self.gon_params = [p for p in self.gon_params if p.num_free() > 0]

    def __call__(self):
        """Perform checks and parameter reduction according to the selected option.

        Returns:
            None
        """

        # In the scan-varying case we can't calculate dB_dp before composing the
        # model, so revert to the original function
        if self._scan_varying:
            self._unit_cell_surplus_reflections = self._surplus_reflections

        # As a special case for detector metrology, try reducing the number of
        # detector parameters if there are too few for some panel group. If this is
        # unsuccessful, fail outright.
        if self._options.detector_reduce:
            self.detector_reduce()

        if self._options.action == "fail":
            self.check_and_fail()

        elif self._options.action == "fix":
            self.check_and_fix()

        elif self._options.action == "remove":
            self.check_and_remove()
