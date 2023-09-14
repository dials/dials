# A port of xfel.mono_simulation.max_like.minimizer
from __future__ import annotations

import logging
import math
import sys
from math import exp, log, pow

from scitbx.array_family import flex
from scitbx.lbfgs import (
    core_parameters,
    exception_handling_parameters,
    run,
    termination_parameters,
)

logger = logging.getLogger(__name__)


class minimizer:
    def __init__(
        self, d_i: flex.double, psi_i: flex.double, eta_rad: float, Deff: float
    ):
        self.safelog = -1.0 + math.log(sys.float_info.max)
        self.S = {"d_i": d_i, "psi_i": psi_i, "eta_rad": eta_rad, "Deff": Deff}
        assert len(d_i) == len(psi_i)
        self.d_i = d_i
        self.psi_i = psi_i
        self.Nobs = len(d_i)
        self.escalate = 10.0  # 10 is a soft switch; 50-100 a hard switch
        self.x = flex.double([log(2.0 / Deff), log(eta_rad)])  # parameters alpha, eta
        self.minimizer = run(
            target_evaluator=self,
            core_params=core_parameters(
                gtol=0.1
                # increasing the accuracy of the line search technique (default=0.9)
                # as suggested by source code.  Otherwise Deff is set unreasonably high
                # and the exponential blows up.
            ),
            termination_params=termination_parameters(
                traditional_convergence_test=False,
                drop_convergence_test_max_drop_eps=1.0e-5,
                min_iterations=0,
                max_iterations=100,
                max_calls=200,
            ),
            exception_handling_params=exception_handling_parameters(
                ignore_line_search_failed_rounding_errors=True,
                ignore_line_search_failed_step_at_lower_bound=True,  # the only change from default
                ignore_line_search_failed_step_at_upper_bound=False,
                ignore_line_search_failed_maxfev=False,
                ignore_line_search_failed_xtol=False,
                ignore_search_direction_not_descent=False,
            ),
        )

        self.x = flex.exp(self.x)

    def functional_only(self, alpha, eta):

        f = 0.0
        for d, psi in zip(self.d_i, self.psi_i):
            psi_model = (d * alpha + eta) / 2.0
            B = self.escalate / psi_model
            expBarg = B * (psi + psi_model)
            expBnegarg = -B * (psi - psi_model)

            if abs(expBarg) > self.safelog or abs(expBnegarg) > self.safelog:
                logger.info(f"XXescalate {self.escalate}")
                logger.info(f"XXpsi_model {psi_model}")
                logger.info(f"XXexp {B} {expBarg}")
                logger.info(f"XXeta {eta}")
                logger.info(f"XXDeff { 1./alpha}")
                logger.info(self.S)
                raise ValueError(
                    f"max likelihood exp argument outside of math domain {expBarg} {expBnegarg}"
                )

            fx = (0.5 / psi_model) / (1 + exp(expBarg)) * (1 + exp(self.escalate))
            gx = 1.0 / (1 + exp(expBnegarg)) * (1 + exp(self.escalate))
            prob = fx * gx
            f -= math.log(prob)
        return f

    def compute_functional_and_gradients(self):
        """The compute_functional_and_gradients() function

        @return Two-tuple of the value of the functional, and an
                <code>n</code>-long vector with the values of the
                gradients at the current position
        """
        if self.x[0] > self.safelog or self.x[1] > self.safelog:
            raise ValueError(
                f"max likelihood current parameters outside of math domain {self.x[0]} {self.x[1]}"
            )
        alpha = exp(self.x[0])
        eta = exp(self.x[1])

        partf_partP0 = 0.0
        partf_partP1 = 0.0

        f = self.functional_only(alpha, eta)

        for d, psi_i in zip(self.d_i, self.psi_i):
            psi_model = (d * alpha + eta) / 2.0
            part_psi_model_partP0 = 0.5 * d * alpha
            part_psi_model_partP1 = 0.5 * eta

            B = self.escalate / psi_model

            if psi_model > 1e100 or part_psi_model_partP0 > 1e100:
                from libtbx.utils import Sorry

                raise Sorry("Model has diverged, cannot continue")

            expB = exp(B * (psi_i + psi_model))
            expBneg = exp(-B * (psi_i - psi_model))

            if expB > 1e100 or expBneg < -1e100 or expBneg > 1e100:
                from libtbx.utils import Sorry

                raise Sorry("Model has diverged, cannot continue")

            Spos = 1.0 + expB
            Sneg = 1.0 + expBneg
            expnu = 1.0 + exp(self.escalate)

            partB_partP0 = (
                -self.escalate / (psi_model * psi_model)
            ) * part_psi_model_partP0
            partB_partP1 = (
                -self.escalate / (psi_model * psi_model)
            ) * part_psi_model_partP1

            partSpos_partP0 = expB * (
                (psi_i + psi_model) * partB_partP0 + B * part_psi_model_partP0
            )
            partSpos_partP1 = expB * (
                (psi_i + psi_model) * partB_partP1 + B * part_psi_model_partP1
            )

            partSneg_partP0 = expBneg * (
                (-psi_i + psi_model) * partB_partP0 + B * part_psi_model_partP0
            )
            partSneg_partP1 = expBneg * (
                (-psi_i + psi_model) * partB_partP1 + B * part_psi_model_partP1
            )

            partG_partP0 = -expnu * pow(Sneg, -2) * partSneg_partP0
            partG_partP1 = -expnu * pow(Sneg, -2) * partSneg_partP1

            Sfac = 2.0 * psi_model * Spos
            partF_partP0 = (
                -expnu
                * pow(Sfac, -2)
                * 2
                * (psi_model * partSpos_partP0 + Spos * part_psi_model_partP0)
            )
            partF_partP1 = (
                -expnu
                * pow(Sfac, -2)
                * 2
                * (psi_model * partSpos_partP1 + Spos * part_psi_model_partP1)
            )

            fx = (0.5 / psi_model) / (Spos) * expnu
            gx = (1.0 / Sneg) * expnu
            prob = fx * gx
            part_prob_partP0 = fx * partG_partP0 + gx * partF_partP0
            part_prob_partP1 = fx * partG_partP1 + gx * partF_partP1

            partf_partP0 -= (1.0 / prob) * part_prob_partP0
            partf_partP1 -= (1.0 / prob) * part_prob_partP1

        return (f, flex.double([partf_partP0, partf_partP1]))

    def fd_compute_functional_and_gradients(self, epsilon=1e-6):
        """The compute_functional_and_gradients() function

        @return Two-tuple of the value of the functional, and an
                <code>n</code>-long vector with the values of the
                gradients at the current position
        """

        alpha = exp(self.x[0])
        eta = exp(self.x[1])
        aplus = exp(self.x[0] + epsilon)
        aminu = exp(self.x[0] - epsilon)
        eplus = exp(self.x[1] + epsilon)
        eminu = exp(self.x[1] - epsilon)

        f = self.functional_only(alpha, eta)

        fd_partf_partalpha = (
            self.functional_only(aplus, eta) - self.functional_only(aminu, eta)
        ) / (2.0 * epsilon)

        fd_partf_parteta = (
            self.functional_only(alpha, eplus) - self.functional_only(alpha, eminu)
        ) / (2.0 * epsilon)

        logger.info(f"{f}, [{fd_partf_partalpha},{fd_partf_parteta}], finite diff\n")

        return (f, flex.double([fd_partf_partalpha, fd_partf_parteta]))
