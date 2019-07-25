#!/usr/bin/env cctbx.python

#
#  Copyright (C) (2016) STFC Rutherford Appleton Laboratory, UK.
#
#  Author: David Waterman.
#
#  This code is distributed under the BSD license, a copy of which is
#  included in the root directory of this package.
#
"""Versions of refinement classes for two theta refinement of the unit cell"""

from __future__ import absolute_import, division, print_function
import logging

logger = logging.getLogger(__name__)

from dials.array_family import flex
from scitbx import matrix
from math import sqrt, pi
from libtbx.table_utils import simple_table
from scitbx.math import five_number_summary
from dials.algorithms.refinement.reflection_manager import ReflectionManager
from dials.algorithms.refinement.prediction.managed_predictors import (
    ExperimentsPredictor,
)
from dials.algorithms.refinement.target import Target
from dials.algorithms.refinement.parameterisation.prediction_parameters import (
    PredictionParameterisation,
)

# constants
RAD2DEG = 180.0 / pi
DEG2RAD = pi / 180.0


class ConstantTwoThetaWeightingStrategy(object):
    def calculate_weights(self, reflections):

        reflections["2theta.weights"] = flex.double(len(reflections), 1)
        return reflections


def calc_2theta(reflections, experiments):
    """Calculate and return 2theta angles in radians"""

    twotheta = flex.double(len(reflections), 0.0)
    for iexp, exp in enumerate(experiments):
        isel = (reflections["id"] == iexp).iselection()
        sub_ref = reflections.select(isel)
        s0 = matrix.col(exp.beam.get_s0())
        for ipanel in range(len(exp.detector)):
            sel = sub_ref["panel"] == ipanel
            panel_ref = sub_ref.select(sel)
            x, y, phi = panel_ref["xyzobs.mm.value"].parts()
            s1 = exp.detector[ipanel].get_lab_coord(flex.vec2_double(x, y))
            s1 = s1 / s1.norms() * s0.length()

            sub_isel = isel.select(sel)
            twotheta.set_selected(sub_isel, s1.angle(s0))

    return twotheta


class TwoThetaReflectionManager(ReflectionManager):
    _weighting_strategy = ConstantTwoThetaWeightingStrategy()

    def __init__(self, *args, **kwargs):

        # call base __init__
        super(TwoThetaReflectionManager, self).__init__(*args, **kwargs)

        # set observed 2theta angles
        self._reflections["2theta_obs.rad"] = calc_2theta(
            self._reflections, self._experiments
        )

        # placeholder for calculated 2theta angles
        self._reflections["2theta_cal.rad"] = flex.double(len(self._reflections), 0.0)

        return

    def print_stats_on_matches(self):

        l = self.get_matches()
        nref = len(l)
        if nref == 0:
            logger.warning(
                "Unable to calculate summary statistics for zero observations"
            )
            return

        twotheta_resid = l["2theta_resid"]
        w_2theta = l["2theta.weights"]

        msg = (
            "\nSummary statistics for {} observations".format(nref)
            + " matched to predictions:"
        )
        header = ["", "Min", "Q1", "Med", "Q3", "Max"]
        rows = []
        row_data = five_number_summary(twotheta_resid)
        rows.append(
            ["2theta_c - 2theta_o (deg)"] + ["%.4g" % (e * RAD2DEG) for e in row_data]
        )
        row_data = five_number_summary(w_2theta)
        rows.append(
            ["2theta weights"] + ["%.4g" % (e * DEG2RAD ** 2) for e in row_data]
        )
        st = simple_table(rows, header)
        logger.info(msg)
        logger.info(st.format())
        logger.info("")


class TwoThetaExperimentsPredictor(ExperimentsPredictor):
    def _predict_one_experiment(self, experiment, reflections):

        B = flex.mat3_double(len(reflections), experiment.crystal.get_B())
        r0 = B * reflections["miller_index"].as_vec3_double()
        r0len = r0.norms()
        wl = experiment.beam.get_wavelength()

        # 2theta = 2 * arcsin( |r0| / (2 * |s0| ) )
        reflections["2theta_cal.rad"] = 2.0 * flex.asin(0.5 * r0len * wl)

        reflections.set_flags(
            flex.size_t(len(reflections)), reflections.flags.predicted
        )


class TwoThetaTarget(Target):
    _grad_names = ["d2theta_dp"]
    rmsd_names = ["RMSD_2theta"]
    rmsd_units = ["rad"]

    def __init__(
        self, experiments, predictor, reflection_manager, prediction_parameterisation
    ):
        Target.__init__(
            self,
            experiments,
            predictor,
            reflection_manager,
            prediction_parameterisation,
        )

        # set the single cutoff for 2theta residual to essentially zero
        self._binsize_cutoffs = [1.0e-6]

        # predict reflections and finalise reflection manager
        self.predict()
        self._reflection_manager.finalise()

        return

    def predict(self):
        """perform reflection prediction for the working reflections and update the
        reflection manager"""

        # get the matches
        reflections = self._reflection_manager.get_obs()

        # reset the 'use' flag for all observations
        self._reflection_manager.reset_accepted_reflections()

        # set twotheta in place
        self._reflection_predictor(reflections)

        # calculate  residuals
        reflections["2theta_resid"] = (
            reflections["2theta_cal.rad"] - reflections["2theta_obs.rad"]
        )
        reflections["2theta_resid2"] = reflections["2theta_resid"] ** 2

        # set used_in_refinement flag to all those that had predictions
        mask = reflections.get_flags(reflections.flags.predicted)
        reflections.set_flags(mask, reflections.flags.used_in_refinement)

        # collect the matches
        self.update_matches(force=True)

        return

    @staticmethod
    def _extract_residuals_and_weights(matches):

        # return residuals and weights as 1d flex.double vectors
        residuals = matches["2theta_resid"]

        weights = matches["2theta.weights"]

        return residuals, weights

    @staticmethod
    def _extract_squared_residuals(matches):

        residuals2 = matches["2theta_resid2"]

        return residuals2

    def _rmsds_core(self, reflections):
        """calculate unweighted RMSDs for the specified reflections"""

        resid_2theta = flex.sum(reflections["2theta_resid2"])
        n = len(reflections)

        rmsds = (sqrt(resid_2theta / n),)
        return rmsds

    def achieved(self):
        """RMSD criterion for target achieved """
        r = self._rmsds if self._rmsds else self.rmsds()

        # reset cached rmsds to avoid getting out of step
        self._rmsds = None

        if r[0] < self._binsize_cutoffs[0]:
            return True
        return False


class TwoThetaPredictionParameterisation(PredictionParameterisation):
    _grad_names = ("d2theta_dp",)

    def __init__(self, *args, **kwargs):
        super(TwoThetaPredictionParameterisation, self).__init__(*args, **kwargs)
        # check that only the unit cell is parameterised
        assert not self._detector_parameterisations
        assert not self._beam_parameterisations
        assert not self._xl_orientation_parameterisations
        assert not self._goniometer_parameterisations
        return

    def _local_setup(self, reflections):

        # we want the wavelength
        self._wavelength = 1.0 / self._s0.norms()

        return

    def _xl_unit_cell_derivatives(self, isel, parameterisation=None, reflections=None):

        # Get required data
        h = self._h.select(isel)
        B = self._B.select(isel)
        wl = self._wavelength.select(isel)

        # get derivatives of the B matrix wrt the parameters
        dB_dxluc_p = [
            None if der is None else flex.mat3_double(len(isel), der.elems)
            for der in parameterisation.get_ds_dp(use_none_as_null=True)
        ]

        d2theta_dp = []

        # loop through the parameters
        for der in dB_dxluc_p:

            if der is None:
                d2theta_dp.append(None)
                continue

            r0 = B * h
            dr0 = der * h
            r0len = r0.norms()
            dr0len = dr0.dot(r0) / r0len

            # 2theta = 2 * arcsin( |r0| / (2 * |s0| ) )
            sintheta = 0.5 * r0len * wl
            fac = 1.0 / flex.sqrt(flex.double(len(wl), 1.0) - sintheta ** 2)
            val = fac * wl * dr0len

            d2theta_dp.append(val)

        return d2theta_dp

    def _grads_xl_unit_cell_loop(self, reflections, results, callback=None):
        """Loop over all crystal unit cell parameterisations, calculate gradients
        and extend the results"""

        # loop over the crystal unit cell parameterisations
        for xlucp in self._xl_unit_cell_parameterisations:

            # Determine (sub)set of reflections affected by this parameterisation
            isel = flex.size_t()
            for exp_id in xlucp.get_experiment_ids():
                isel.extend(self._experiment_to_idx[exp_id])

            # Extend derivative vectors for this crystal unit cell parameterisation
            results = self._extend_gradient_vectors(
                results, self._nref, xlucp.num_free(), keys=self._grad_names
            )

            if len(isel) == 0:
                # if no reflections are in this experiment, skip calculation of
                # gradients, but must still process null gradients by a callback
                if callback is not None:
                    for iparam in range(xlucp.num_free()):
                        results[self._iparam] = callback(results[self._iparam])
                        self._iparam += 1
                else:
                    self._iparam += xlucp.num_free()
                continue

            d2theta_dp = self._xl_unit_cell_derivatives(
                isel, parameterisation=xlucp, reflections=reflections
            )

            for d2theta in d2theta_dp:
                if d2theta is not None:
                    results[self._iparam][self._grad_names[0]].set_selected(
                        isel, d2theta
                    )

                # increment the parameter index pointer
                self._iparam += 1

        return results
