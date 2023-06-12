from __future__ import annotations

import logging
import random
import textwrap
from math import log, pi, sqrt
from typing import List

import numpy as np
from numpy.linalg import det, inv, norm

from dxtbx import flumpy
from scitbx import linalg, matrix

from dials.algorithms.profile_model.ellipsoid import mosaicity_from_eigen_decomposition
from dials.algorithms.profile_model.ellipsoid.model import (
    compute_change_of_basis_operation,
)
from dials.algorithms.profile_model.ellipsoid.parameterisation import (
    ReflectionModelState,
)
from dials.array_family import flex
from dials.util import tabulate
from dials_algorithms_profile_model_ellipsoid_ext import reflection_statistics, rse

logger = logging.getLogger("dials")

flex.set_random_seed(0)
random.seed(0)


class BadSpotForIntegrationException(Exception):
    pass


def compute_dSbar(S: np.array, dS: np.array) -> np.array:
    # dS & S are 3x3 arrays. Returns a 2x2 array
    S12 = S[0:2, 2].reshape(2, 1)
    S21 = S[2, 0:2].reshape(1, 2)
    S22 = S[2, 2]

    dS11 = dS[0:2, 0:2]
    dS12 = dS[0:2, 2].reshape(2, 1)
    dS21 = dS[2, 0:2].reshape(1, 2)
    dS22 = dS[2, 2]

    S22_inv = 1 / S22

    A = dS11
    B = np.matmul(S12 * S22_inv * dS22 * S22_inv, S21)
    C = np.matmul(S12 * S22_inv, dS21)
    D = np.matmul(dS12 * S22_inv, S21)
    return A + B - (C + D)


def compute_dmbar(S: np.array, dS: np.array, dmu: np.array, epsilon: float) -> np.array:
    # S, dS are is a 3x3 array, dmu is 3x1 array
    dmu = dmu.reshape(3, 1)  # 3x1 array

    S12 = S[0:2, 2].reshape(2, 1)
    S22 = S[2, 2]

    dS12 = dS[0:2, 2].reshape(2, 1)
    dS22 = dS[2, 2]

    S22_inv = 1 / S22

    dmu1 = dmu[0:2, 0].reshape(2, 1)
    dmu2 = dmu[2, 0]
    dep = -dmu2

    A = dmu1
    B = dS12 * S22_inv * epsilon
    C = -S12 * S22_inv * dS22 * S22_inv * epsilon
    D = S12 * S22_inv * dep
    return A + B + C + D


class ConditionalDistribution(object):
    """
    A class to compute useful stuff about the conditional distribution

    """

    def __init__(self, norm_s0, mu, dmu, S, dS):
        # norm_s0 is a float i.e. norm(s0)
        self._mu = mu  # 3x1 array
        self._dmu = dmu  # 3 x n array
        self._S = S  # 3x3 array
        self._dS = dS  # 3 x 3 x n array

        # Partition the covariance matrix
        S11 = S[0:2, 0:2]
        S12 = S[0:2, 2].reshape(2, 1)
        S21 = S[2, 0:2].reshape(1, 2)
        S22 = S[2, 2]

        # The partitioned mean vector
        mu1 = mu[0:2, 0].reshape(2, 1)
        mu2 = mu[2, 0]
        # a = norm(s0)

        # The epsilon
        self.epsilon = norm_s0 - mu2

        # Compute the conditional mean
        self._mubar = mu1 + S12 * (1 / S22) * self.epsilon
        assert self._mubar.shape == (2, 1)

        # Compute the conditional covariance matrix
        self._Sbar = S11 - np.matmul(S12 * (1 / S22), S21)
        assert self._Sbar.shape == (2, 2)

        # Set to None and compute on demand
        self.dSbar = None
        self.dmbar = None
        self.d2Sbar = None
        self.d2mbar = None

    def mean(self) -> np.array:
        """
        Return the conditional mean (a 2x1 array)

        """
        return self._mubar

    def sigma(self) -> np.array:
        """
        Return the conditional sigma (a 2x2 array)

        """
        return self._Sbar

    def first_derivatives_of_sigma(self) -> List[np.array]:
        """
        Return the marginal first derivatives (as a list of 2x2 arrays)

        """
        if self.dSbar is None:
            self.dSbar = [
                compute_dSbar(self._S, self._dS[:, :, i])
                for i in range(self._dS.shape[2])
            ]

        return self.dSbar

    def first_derivatives_of_mean(self) -> List[np.array]:
        """
        Return the marginal first derivatives (a list of 2x1 arrays)

        """
        if self.dmbar is None:
            self.dmbar = [
                compute_dmbar(self._S, self._dS[:, :, i], self._dmu[:, i], self.epsilon)
                for i in range(self._dS.shape[2])
            ]

        return self.dmbar


def rotate_vec3_double(R, A):
    """
    Helper function to rotate an array of matrices

    """
    return np.einsum("ij,jk->ik", R, A)


def rotate_mat3_double(R, A):
    """
    Helper function to rotate an array of matrices

    """
    return np.einsum("ij,jkv,kl->ilv", R, A, R.T)


class ReflectionLikelihood(object):
    def __init__(self, model, s0, sp, h, ctot, mobs, sobs, panel_id=0):

        # Save stuff
        modelstate = ReflectionModelState(model, s0, h)
        self.modelstate = modelstate
        self.s0 = s0.reshape(3, 1)
        self.norm_s0 = norm(s0)
        self.sp = sp.reshape(3, 1)
        self.h = np.array([h], dtype=np.float64).reshape(3, 1)
        self.ctot = ctot
        self.mobs = mobs.reshape(2, 1)
        self.sobs = sobs
        self.panel_id = panel_id

        # Compute the change of basis
        self.R = compute_change_of_basis_operation(self.s0, self.sp)  # const
        self.R_cctbx = matrix.sqr(flex.double(self.R.tolist()))
        s2 = self.s0 + self.modelstate.get_r()
        # Rotate the mean vector
        self.mu = np.matmul(self.R, s2)
        self.S = np.matmul(
            np.matmul(self.R, modelstate.mosaicity_covariance_matrix), self.R.T
        )  # const when not refining mosaicity
        self.dS = rotate_mat3_double(
            self.R, modelstate.get_dS_dp()
        )  # const when not refining mosaicity?
        self.dmu = rotate_vec3_double(
            self.R, modelstate.get_dr_dp()
        )  # const when not refining uc/orientation?
        # Construct the conditional distribution
        self.conditional = ConditionalDistribution(
            self.norm_s0, self.mu, self.dmu, self.S, self.dS
        )

    def update(self):

        # The s2 vector
        s2 = self.s0 + self.modelstate.get_r()
        # Rotate the mean vector
        self.mu = np.matmul(self.R, s2)

        # Rotate the covariance matrix
        if not self.modelstate.state.is_mosaic_spread_fixed:
            self.S = np.matmul(
                np.matmul(self.R, self.modelstate.mosaicity_covariance_matrix), self.R.T
            )  # const when not refining mosaicity

        # Rotate the first derivative matrices
        if not self.modelstate.state.is_mosaic_spread_fixed:
            self.dS = rotate_mat3_double(
                self.R, self.modelstate.get_dS_dp()
            )  # const when not refining mosaicity?

        # Rotate the first derivative of s2
        if (not self.modelstate.state.is_unit_cell_fixed) or not (
            self.modelstate.state.is_orientation_fixed
        ):
            self.dmu = rotate_vec3_double(
                self.R, self.modelstate.get_dr_dp()
            )  # const when not refining uc/orientation?

        # Construct the conditional distribution
        self.conditional = ConditionalDistribution(
            self.norm_s0, self.mu, self.dmu, self.S, self.dS
        )

    def log_likelihood(self):
        """
        Compute the log likelihood for the reflection

        """

        # Get data
        ctot = self.ctot
        mobs = self.mobs
        Sobs = self.sobs

        # Get info about the marginal
        S22 = self.S[2, 2]
        S22_inv = 1 / S22
        mu2 = self.mu[2, 0]

        # Get info about the conditional
        Sbar = self.conditional.sigma()
        mubar = self.conditional.mean()
        Sbar_inv = inv(Sbar)
        Sbar_det = det(Sbar)

        # Weights for marginal and conditional components
        m_w = ctot
        c_w = ctot

        # Compute the marginal likelihood
        m_d = self.norm_s0 - mu2
        m_lnL = m_w * (log(S22) + S22_inv * m_d**2)

        # Compute the conditional likelihood
        c_d = mobs - mubar
        c_lnL = c_w * (
            log(Sbar_det)
            + np.trace(np.matmul(Sbar_inv, (Sobs + np.matmul(c_d, c_d.T))))
        )

        # Return the joint likelihood
        jLL = -0.5 * (m_lnL + c_lnL)
        return jLL

    def first_derivatives(self):
        """
        Compute the first derivatives

        """
        # Get data
        ctot = self.ctot
        mobs = self.mobs  # 2x1 array
        Sobs = self.sobs

        # Get info about marginal distribution
        S22 = self.S[2, 2]
        dS22_vec = self.dS[2, 2, :]  # for i in range(self.dS.shape[2])]
        S22_inv = 1 / S22
        mu2 = self.mu[2, 0]

        # Get info about conditional distribution
        Sbar = self.conditional.sigma()  # 3x3 array
        mubar = self.conditional.mean()  # 2x1 array
        dSbar = self.conditional.first_derivatives_of_sigma()  # list of 2x2 arrays
        dmbar = self.conditional.first_derivatives_of_mean()  # list of 2x1 arrays
        Sbar_inv = inv(Sbar)

        # The distance from the ewald sphere
        epsilon = self.norm_s0 - mu2
        c_d = mobs - mubar  # 2x1 array

        # Weights for marginal and conditional components
        m_w = ctot
        c_w = ctot

        # Compute the derivative wrt parameter i
        I = np.array([[1.0, 0], [0, 1.0]], dtype=np.float64).reshape(2, 2)

        V1 = Sobs + np.matmul(c_d, c_d.T)
        V2 = I - np.matmul(Sbar_inv, V1)

        dSbar = np.array(dSbar)
        dmbar_vec = np.array(dmbar)

        V_vec = np.einsum("ij,ljk->ikl", Sbar_inv, dSbar)
        V_vec = c_w * np.einsum("ijl,ji->l", V_vec, V2)
        dep = -self.dmu[2, :]
        U_vec = m_w * (
            S22_inv * dS22_vec * (1.0 - S22_inv * epsilon**2)
            + 2 * S22_inv * epsilon * dep
        )
        W_vec = np.einsum("ij,lkj->ikl", c_d, dmbar_vec)
        W_vec = -2.0 * c_w * np.einsum("ij,jil->l", Sbar_inv, W_vec)

        dL = -0.5 * (U_vec + V_vec + W_vec)
        return dL

    def fisher_information(self):
        """
        Compute the fisher information

        """
        ctot = self.ctot

        # Get info about marginal distribution
        S22 = self.S[2, 2]
        dS22 = [self.dS[2, 2, i] for i in range(self.dS.shape[2])]
        S22_inv = 1 / S22

        # Get info about conditional distribution
        Sbar = self.conditional.sigma()
        dSbar = self.conditional.first_derivatives_of_sigma()  # list of 2x2 arrays
        dmbar = self.conditional.first_derivatives_of_mean()  # list of 2x1 arrays
        Sbar_inv = inv(Sbar)
        dmu = self.dmu

        # Weights for marginal and conditional components
        m_w = ctot
        c_w = ctot

        # Compute the fisher information wrt parameter i j
        I = flex.double(flex.grid(len(dS22), len(dS22)))

        for j in range(len(dS22)):
            for i in range(len(dS22)):
                U = S22_inv * dS22[j] * S22_inv * dS22[i]
                V = np.trace(
                    np.matmul(
                        np.matmul(
                            np.matmul(
                                Sbar_inv,
                                dSbar[j],
                            ),
                            Sbar_inv,
                        ),
                        dSbar[i],
                    )
                )
                W = 2 * np.trace(
                    np.matmul(
                        np.matmul(
                            Sbar_inv,
                            dmbar[i],
                        ),
                        dmbar[j].T,
                    )
                )
                X = 2 * dmu[2, i] * S22_inv * dmu[2, j]
                I[j, i] = 0.5 * c_w * (V + W) + 0.5 * m_w * (U + X)

        return I


class MaximumLikelihoodTarget(object):
    def __init__(
        self, model, s0, sp_list, h_list, ctot_list, mobs_list, sobs_list, panel_ids
    ):

        # Check input
        assert len(h_list) == sp_list.shape[-1]
        assert len(h_list) == ctot_list.shape[-1]
        assert len(h_list) == mobs_list.shape[-1]
        assert len(h_list) == sobs_list.shape[-1]

        # Save the model
        self.model = model

        # Compute the change of basis for each reflection
        self.data = []
        for i in range(len(h_list)):
            self.data.append(
                ReflectionLikelihood(
                    model,
                    s0,
                    sp_list[:, i],
                    matrix.col(h_list[i]),
                    ctot_list[i],
                    mobs_list[:, i],
                    sobs_list[:, :, i],
                    panel_ids[i],
                )
            )

    def update(self):
        for d in self.data:
            d.modelstate.update()  # update the ReflectionModelState
            d.update()  # update the ReflectionLikelihood

    def mse(self):
        """
        The MSE in local reflection coordinates

        """
        mse = 0
        for i in range(len(self.data)):
            mbar = self.data[i].conditional.mean()
            xobs = self.data[i].mobs
            mse += np.dot((xobs - mbar).T, xobs - mbar)
        mse /= len(self.data)
        return mse

    def rmsd(self):
        """
        The RMSD in pixels

        """
        mse_x = 0.0
        mse_y = 0.0
        for i in range(len(self.data)):
            R = self.data[i].R_cctbx
            mbar = tuple(self.data[i].conditional.mean().flatten())
            xobs = tuple(self.data[i].mobs.flatten())
            norm_s0 = self.data[i].norm_s0
            rse_i = rse(R, mbar, xobs, norm_s0, self.model.experiment.detector)
            mse_x += rse_i[0]
            mse_y += rse_i[1]
        mse_x /= len(self.data)
        mse_y /= len(self.data)
        return np.sqrt(np.array([mse_x, mse_y]))

    def log_likelihood(self):
        """
        The joint log likelihood

        """
        return sum(d.log_likelihood() for d in self.data)

    def jacobian(self):
        """
        Return the Jacobian

        """
        return flex.double([list(d.first_derivatives()) for d in self.data])

    def first_derivatives(self):
        """
        The joint first derivatives

        """
        dL = 0
        for d in self.data:
            dL += d.first_derivatives()
        return dL

    def fisher_information(self):
        """
        The joint fisher information

        """
        return sum(d.fisher_information() for d in self.data)


def line_search(func, x, p, tau=0.5, delta=1.0, tolerance=1e-7):
    """
    Perform a line search
    :param func The function to minimize
    :param x The initial position
    :param p The direction to search
    :param tau: The backtracking parameter
    :param delta: The initial step
    :param tolerance: The algorithm tolerance
    :return: The amount to move in the given direction

    """
    fa = func(x)
    if p.length() < 1:
        min_delta = tolerance
    else:
        min_delta = tolerance / p.length()

    while delta > min_delta:
        try:
            fb = func(x + delta * p)
            if fb <= fa:
                return delta
        except Exception:
            pass
        delta *= tau
    return 0


def gradient_descent(f, df, x0, max_iter=1000, tolerance=1e-10):
    """
    Find the minimum using gradient descent and a line search
    :param f The function to minimize
    :param df The function to compute derivatives
    :param x0 The initial position
    :param max_iter: The maximum number of iterations
    :param tolerance: The algorithm tolerance
    :return: The amount to move in the given direction

    """
    delta = 0.5
    for it in range(max_iter):
        p = -matrix.col(df(x0))
        delta = line_search(f, x0, p, delta=min(1.0, delta * 2), tolerance=tolerance)
        x = x0 + delta * p
        assert f(x) <= f(x0)
        if (x - x0).length() < tolerance:
            break
        x0 = x
    return x


class FisherScoringMaximumLikelihoodBase(object):
    """
    A class to solve maximum likelihood equations using fisher scoring

    """

    def __init__(self, x0, max_iter=1000, tolerance=1e-7, LL_tolerance=1e-6):
        """
        Configure the algorithm

        :param x0: The initial parameter estimates
        :param max_iter: The maximum number of iterations
        :param tolerance: The parameter tolerance

        """
        self.x0 = matrix.col(x0)
        self.max_iter = max_iter
        self.tolerance = tolerance
        self.LL_tolerance = LL_tolerance

    def solve(self):
        """
        Find the maximum likelihood estimate

        """
        x0 = self.x0

        # Loop through the maximum number of iterations
        for it in range(self.max_iter):
            # Compute the derivative and fisher information at x0
            S, I = self.score_and_fisher_information(x0)

            # Solve the update equation to get direction
            p = matrix.col(self.solve_update_equation(S, I))

            # Perform a line search to ensure that each step results in an increase the
            # in log likelihood. In the rare case where the update does not result in an
            # increase in the likelihood (only observed for absurdly small samples
            # e.g. 2 reflections or when 1 parameter approaches zero) do an iteration
            # of gradient descent
            delta = self.line_search(x0, p)
            if delta > 0:
                x = x0 + delta * p
            else:
                x = self.gradient_search(x0)

            # Call an update
            self.callback(x)
            # Break the loop if the parameters change less than the tolerance
            if (x - x0).length() < self.tolerance:
                break
            if self.test_LL_convergence():
                break

            # Update the parameter
            x0 = x

        # Save the parameters
        self.num_iter = it + 1
        self.parameters = x

    def test_LL_convergence(self):
        try:
            l1 = self.history[-1]["likelihood"]
            l2 = self.history[-2]["likelihood"]
        except IndexError:
            return False

        test = abs(l1 - l2) < self.LL_tolerance
        return test

    def solve_update_equation(self, S, I):
        """
        Solve the update equation using cholesky decomposition
        :param S: The score
        :param I: The fisher information
        :return: The parameter delta

        """

        # Construct triangular matrix
        LL = flex.double()
        for j in range(len(S)):
            for i in range(j + 1):
                LL.append(I[j * len(S) + i])

        # Perform the decomposition
        ll = linalg.l_l_transpose_cholesky_decomposition_in_place(LL)
        p = flex.double(S)
        return ll.solve(p)

    def line_search(self, x, p, tau=0.5, delta=1.0, tolerance=1e-7):
        """
        Perform a line search
        :param x The initial position
        :param p The direction to search
        :return: The amount to move in the given direction

        """

        def f(x):
            return -self.log_likelihood(x)

        return line_search(f, x, p, tolerance=self.tolerance)

    def gradient_search(self, x0):
        """
        Find the minimum using gradient descent and a line search
        :param x0 The initial position
        :return: The amount to move in the given direction

        """

        def f(x):
            return -self.log_likelihood(x)

        def df(x):
            return -self.score(x)

        return gradient_descent(f, df, x0, max_iter=1, tolerance=self.tolerance)


class FisherScoringMaximumLikelihood(FisherScoringMaximumLikelihoodBase):
    """
    A class to solve the maximum likelihood equations

    """

    def __init__(
        self,
        model,
        s0,
        sp_list,
        h_list,
        ctot_list,
        mobs_list,
        sobs_list,
        panel_ids,
        max_iter=1000,
        tolerance=1e-7,
        LL_tolerance=1e-6,
    ):
        """
        Initialise the algorithm:

        """
        # Initialise the super class
        super(FisherScoringMaximumLikelihood, self).__init__(
            model.active_parameters,
            max_iter=max_iter,
            tolerance=tolerance,
            LL_tolerance=LL_tolerance,
        )

        # Save the parameterisation
        self.model = model

        # Save some stuff
        self.s0 = s0
        self.sp_list = sp_list
        self.h_list = h_list
        self.ctot_list = ctot_list
        self.mobs_list = mobs_list
        self.sobs_list = sobs_list
        self.panel_ids = panel_ids

        # Store the parameter history
        self.history = []

        self._ml_target = MaximumLikelihoodTarget(
            self.model,
            self.s0,
            self.sp_list,
            self.h_list,
            self.ctot_list,
            self.mobs_list,
            self.sobs_list,
            self.panel_ids,
        )

        # Print initial
        self.callback(self.model.active_parameters)

    def log_likelihood(self, x):
        """
        :param x: The parameter estimate
        :return: The log likelihood at x

        """
        self.model.active_parameters = x
        self._ml_target.update()
        return self._ml_target.log_likelihood()

    def score(self, x):
        """
        :param x: The parameter estimate
        :return: The score at x

        """
        self.model.active_parameters = x
        self._ml_target.update()
        return flumpy.from_numpy(self._ml_target.first_derivatives())

    def score_and_fisher_information(self, x):
        """
        :param x: The parameter estimate
        :return: The score and fisher information at x

        """
        self.model.active_parameters = x
        self._ml_target.update()
        S = flumpy.from_numpy(self._ml_target.first_derivatives())
        I = self._ml_target.fisher_information()
        return S, I

    def mse(self, x):
        """
        :param x: The parameter estimate
        :return: The MSE at x

        """
        return self._ml_target.mse()

    def rmsd(self, x):
        """
        :param x: The parameter estimate
        :return: The RMSD at x

        """
        return self._ml_target.rmsd()

    def jacobian(self, x):
        """
        :param x: The parameter estimate
        :return: The Jacobian at x

        """
        return self._ml_target.jacobian()

    def condition_number(self, x):
        """
        The condition number of the Jacobian

        """
        from scitbx.linalg.svd import real as svd_real

        svd = svd_real(self.jacobian(x), False, False)
        return max(svd.sigma) / min(svd.sigma)

    def correlation(self, x):
        """
        The correlation of the Jacobian

        """
        J = self.jacobian(x)
        C = flex.double(flex.grid(J.all()[1], J.all()[1]))
        for j in range(C.all()[0]):
            for i in range(C.all()[1]):
                a = J[:, i : i + 1].as_1d()
                b = J[:, j : j + 1].as_1d()
                C[j, i] = flex.linear_correlation(a, b).coefficient()
        return C

    def callback(self, x):
        """
        Handle and update in parameter values

        """
        self.model.active_parameters = x
        self._ml_target.update()
        lnL = self._ml_target.log_likelihood()
        mse = self._ml_target.mse()
        rmsd = self._ml_target.rmsd()

        # Get the unit cell
        unit_cell = self.model.unit_cell.parameters()

        # Get some matrices
        U = self.model.U_matrix.flatten()
        M = (
            self.model._M_parameterisation.sigma().flatten()
        )  # mosaicity_covariance_matrix.flatten()

        # Print some information
        format_string1 = "  Unit cell: (%.3f, %.3f, %.3f, %.3f, %.3f, %.3f)"
        format_string2 = "  | % .6f % .6f % .6f |"
        format_string3 = "  | % .2e % .2e % .2e |"

        logger.info(f"\nIteration: {len(self.history)}")
        if not self.model.is_unit_cell_fixed:
            logger.info("\n" + format_string1 % unit_cell)
        if not self.model.is_orientation_fixed:
            logger.info(
                "\n".join(
                    [
                        "",
                        "  U matrix (orientation)",
                        format_string2 % tuple(U[0:3]),
                        format_string2 % tuple(U[3:6]),
                        format_string2 % tuple(U[6:9]),
                    ]
                )
            )
        if not self.model.is_mosaic_spread_fixed:
            logger.info(
                "\n".join(
                    [
                        "",
                        "  Sigma M",
                        format_string3 % tuple(M[0:3]),
                        format_string3 % tuple(M[3:6]),
                        format_string3 % tuple(M[6:9]),
                    ]
                )
            )
            if self.model.is_mosaic_spread_angular:
                MA = self.model._M_parameterisation.sigma_A().flatten()
                logger.info(
                    "\n".join(
                        [
                            "",
                            "  Sigma M",
                            format_string3 % tuple(MA[0:3]),
                            format_string3 % tuple(MA[3:6]),
                            format_string3 % tuple(MA[6:9]),
                        ]
                    )
                )

        logger.info(
            "\n".join(
                [
                    "",
                    "  ln(L) = %f" % lnL,
                    "",
                    "  R.M.S.D (local) = %.5g" % sqrt(mse),
                    "",
                    "  R.M.S.D (pixel): X = %.3f, Y = %.3f" % tuple(rmsd),
                ]
            )
        )

        # Append the parameters to the history
        self.history.append(
            {
                "parameters": list(x),
                "likelihood": lnL,
                "unit_cell": unit_cell,
                "orientation": list(U),
                "rlp_mosaicity": list(M),
                "rmsd": tuple(rmsd),
            }
        )


class Refiner(object):
    """
    High level profile refiner class that handles book keeping etc

    """

    def __init__(self, state, data, max_iter=1000, LL_tolerance=1e-6):
        """
        Set the data and initial parameters

        """
        self.s0 = data.s0
        self.h_list = data.h_list
        self.sp_list = data.sp_list
        self.ctot_list = data.ctot_list
        self.mobs_list = data.mobs_list
        self.sobs_list = data.sobs_list
        self.panel_ids = data.panel_ids
        self.state = state
        self.history = []
        self.max_iter = max_iter
        self.LL_tolerance = LL_tolerance

    def refine(self):
        """
        Perform the profile refinement

        """
        self.refine_fisher_scoring()

    def refine_fisher_scoring(self):
        """
        Perform the profile refinement

        """

        # Print information
        logger.info("\nComponents to refine:")
        logger.info(" Orientation:       %s" % (not self.state.is_orientation_fixed))
        logger.info(" Unit cell:         %s" % (not self.state.is_unit_cell_fixed))
        logger.info(" RLP mosaicity:     %s" % (not self.state.is_mosaic_spread_fixed))
        logger.info(
            " Wavelength spread: %s\n" % (not self.state.is_wavelength_spread_fixed)
        )

        # Initialise the algorithm
        self.ml = FisherScoringMaximumLikelihood(
            self.state,
            self.s0,
            self.sp_list,
            self.h_list,
            self.ctot_list,
            self.mobs_list,
            self.sobs_list,
            self.panel_ids,
            max_iter=self.max_iter,
            LL_tolerance=self.LL_tolerance,
        )

        # Solve the maximum likelihood equations
        self.ml.solve()

        # Get the parameters
        self.parameters = flex.double(self.ml.parameters)

        # set the parameters
        self.state.active_parameters = self.parameters

        # Print summary table of refinement.
        rows = []
        headers = ["Iteration", "likelihood", "RMSD (pixel) X,Y"]
        for i, h in enumerate(self.ml.history):
            l = h["likelihood"]
            rmsd = h["rmsd"]
            rows.append([str(i), f"{l:.4f}", f"{rmsd[0]:.3f}, {rmsd[1]:.3f}"])
        logger.info(
            "\nRefinement steps:\n\n" + textwrap.indent(tabulate(rows, headers), " ")
        )

        # Print the eigen values and vectors of sigma_m
        if not self.state.is_mosaic_spread_fixed:
            logger.info("\nDecomposition of Sigma_M:")
            print_eigen_values_and_vectors_static(
                matrix.sqr(
                    flumpy.from_numpy(self.state._M_parameterisation.sigma().flatten())
                )
            )
            if self.state.is_mosaic_spread_angular:
                logger.info(
                    print_eigen_values_and_vectors_angular(
                        matrix.sqr(
                            flumpy.from_numpy(
                                self.state._M_parameterisation.sigma_A()[
                                    :2, :2
                                ].flatten()
                            )
                        )
                    )
                )

        # Save the history
        self.history = self.ml.history

        # Return the optimizer
        return self.ml

    def correlation(self):
        """
        Return the correlation matrix between parameters

        """
        return self.ml.correlation(self.state.active_parameters)

    def labels(self):
        """
        Return parameter labels

        """
        return self.state.parameter_labels


class RefinerData(object):
    """
    A class for holding the data needed for the profile refinement

    """

    def __init__(self, s0, sp_list, h_list, ctot_list, mobs_list, sobs_list, panel_ids):
        """
        Init the data

        ctot_list is a list of total counts per reflection

        """
        self.s0 = s0
        self.sp_list = sp_list
        self.h_list = h_list
        self.ctot_list = ctot_list
        self.mobs_list = mobs_list
        self.sobs_list = sobs_list
        self.panel_ids = panel_ids

    @classmethod
    def from_reflections(self, experiment, reflections):
        """
        Generate the required data from the reflections

        """

        # Get the beam vector
        s0 = np.array([experiment.beam.get_s0()], dtype=np.float64).reshape(3, 1)

        # Get the reciprocal lattice vector
        h_list = reflections["miller_index"]

        # Initialise the list of observed intensities and covariances
        sp_list = np.zeros(shape=(3, len(h_list)))
        ctot_list = np.zeros(shape=(len(h_list)))
        mobs_list = np.zeros(shape=(2, len(h_list)))
        Sobs_list = np.zeros(shape=(2, 2, len(h_list)))

        logger.info(
            "Computing observed covariance for %d reflections" % len(reflections)
        )
        s0_length = norm(s0)

        sbox = reflections["shoebox"]
        for r, (panel_id, xyz) in enumerate(
            zip(reflections["panel"], reflections["xyzobs.px.value"])
        ):
            panel = experiment.detector[panel_id]
            sp, ctot, xbar, Sobs = reflection_statistics(
                panel, xyz, s0_length, experiment.beam.get_s0(), sbox[r]
            )

            # Check we have a sensible number of counts
            if ctot <= 0:
                raise BadSpotForIntegrationException(
                    "Strong spot found with <= 0 counts! Check spotfinding results"
                )

            if (Sobs[0] <= 0) or (Sobs[3] <= 0):
                raise BadSpotForIntegrationException(
                    "Strong spot variance <= 0. Check spotfinding results"
                )

            # Add to the lists
            sp_list[:, r] = sp  # [:, 0]
            ctot_list[r] = ctot
            mobs_list[:, r] = xbar  # [:, 0]
            Sobs_list[:, :, r] = np.array([Sobs[0:2], Sobs[2:]], dtype=np.float64)

        # Print some information
        logger.info("")
        logger.info(
            "I_min = %.2f, I_max = %.2f" % (np.min(ctot_list), np.max(ctot_list))
        )

        # Sometimes a single reflection might have an enormouse intensity for
        # whatever reason and since we weight by intensity, this can cause the
        # refinement to be dominated by these reflections. Therefore, if the
        # intensity is greater than some value, damp the weighting accordingly
        def damp_outlier_intensity_weights(ctot_list):
            n = ctot_list.size
            sorted_ctot = np.sort(ctot_list)
            Q1 = sorted_ctot[n // 4]
            Q2 = sorted_ctot[n // 2]
            Q3 = sorted_ctot[3 * n // 4]
            IQR = Q3 - Q1
            T = Q3 + 1.5 * IQR
            logger.info(f"Median I = {Q2:.2f}\nQ1/Q3 I = {Q1:.2f}, {Q3:.2f}")
            logger.info(f"Damping effect of intensities > {T:.2f}")
            ndamped = 0
            for i, ctot in enumerate(ctot_list):
                if ctot > T:
                    logger.debug(f"Damping {ctot:.2f}")
                    ctot_list[i] = T
                    ndamped += 1
            logger.info(f"Damped {ndamped}/{n} reflections")
            return ctot_list

        ctot_list = damp_outlier_intensity_weights(ctot_list)

        # Print the mean covariance
        Smean = np.mean(Sobs_list, axis=2)
        logger.info("")
        logger.info("Mean observed covariance:")
        logger.info(print_matrix_np(Smean))
        print_eigen_values_and_vectors_of_observed_covariance(Smean, s0)

        # Compute the distance from the Ewald sphere
        epsilon = flex.double(
            s0_length - matrix.col(s).length() for s in reflections["s2"]
        )
        mv = flex.mean_and_variance(epsilon)
        logger.info("")
        logger.info("Mean distance from Ewald sphere: %.3g" % mv.mean())
        logger.info(
            "Variance in distance from Ewald sphere: %.3g"
            % mv.unweighted_sample_variance()
        )

        # Return the profile refiner data
        return RefinerData(
            s0, sp_list, h_list, ctot_list, mobs_list, Sobs_list, reflections["panel"]
        )


def print_eigen_values_and_vectors_of_observed_covariance(A, s0):
    """
    Print the eigen values and vectors of a matrix

    """

    # Compute the eigen decomposition of the covariance matrix
    A = matrix.sqr(flumpy.from_numpy(A))
    s0 = matrix.col(flumpy.from_numpy(s0))
    eigen_decomposition = linalg.eigensystem.real_symmetric(A.as_flex_double_matrix())
    Q = matrix.sqr(eigen_decomposition.vectors())
    L = matrix.diag(eigen_decomposition.values())

    # Print the matrix eigen values
    logger.info(f"\nEigen Values:\n{print_matrix(L, indent=2)}\n")
    logger.info(f"\nEigen Vectors:\n{print_matrix(Q, indent=2)}\n")

    logger.info("Observed covariance in degrees equivalent units")
    logger.info("C1: %.5f degrees" % (sqrt(L[0]) * (180.0 / pi) / s0.length()))
    logger.info("C2: %.5f degrees" % (sqrt(L[3]) * (180.0 / pi) / s0.length()))


def print_eigen_values_and_vectors_static(A):
    """
    Print the eigen values and vectors of a matrix

    """

    # Compute the eigen decomposition of the covariance matrix
    eigen_decomposition = linalg.eigensystem.real_symmetric(A.as_flex_double_matrix())
    eigen_values = eigen_decomposition.values()

    # Print the matrix eigen values
    logger.info(
        f"\n Eigen Values:\n{print_matrix(matrix.diag(eigen_values), indent=2)}\n"
    )
    logger.info(
        f"\n Eigen Vectors:\n{print_matrix(matrix.sqr(eigen_decomposition.vectors()), indent=2)}\n"
    )
    logger.info(
        f"""
 Invariant crystal mosaicity:
 M1 : {eigen_values[0]**0.5:.5f} Å⁻¹
 M2 : {eigen_values[1]**0.5:.5f} Å⁻¹
 M3 : {eigen_values[2]**0.5:.5f} Å⁻¹
"""
    )


def print_eigen_values_and_vectors_angular(A):
    """
    Print the eigen values and vectors of a matrix

    """

    # Compute the eigen decomposition of the covariance matrix
    eigen_decomposition = linalg.eigensystem.real_symmetric(A.as_flex_double_matrix())
    eigen_values = eigen_decomposition.values()

    # Print the matrix eigen values
    logger.info(
        f"\n Eigen Values:\n{print_matrix(matrix.diag(eigen_values), indent=2)}\n"
    )
    logger.info(
        f"\n Eigen Vectors:\n{print_matrix(matrix.sqr(eigen_decomposition.vectors()), indent=2)}\n"
    )

    mosaicity = mosaicity_from_eigen_decomposition(eigen_values)
    logger.info(
        """
 Angular Mosaicity in degrees equivalent units:\n"""
        + "\n".join(f" M{i+1} : {m:.5f} degrees" for i, m in enumerate(mosaicity))
    )


def print_matrix_np(A, fmt="%.3g", indent=0):
    """
    Pretty print matrix

    """
    t = [fmt % a for a in A.flatten()]
    l = [len(tt) for tt in t]
    max_l = max(l)
    fmt = "%" + ("%d" % (max_l + 1)) + "s"
    prefix = " " * indent
    lines = []
    for j in range(A.shape[0]):
        line = ""
        for i in range(A.shape[1]):
            line += fmt % t[i + j * A.shape[1]]
        lines.append("%s|%s|" % (prefix, line))
    return "\n".join(lines)


def print_matrix(A, fmt="%.3g", indent=0):
    """
    Pretty print matrix

    """
    t = [fmt % a for a in A]
    l = [len(tt) for tt in t]
    max_l = max(l)
    fmt = "%" + ("%d" % (max_l + 1)) + "s"
    prefix = " " * indent
    lines = []
    for j in range(A.n[0]):
        line = ""
        for i in range(A.n[1]):
            line += fmt % t[i + j * A.n[1]]
        lines.append("%s|%s|" % (prefix, line))
    return "\n".join(lines)
