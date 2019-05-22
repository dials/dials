from __future__ import division

import logging
from math import log, pi, sqrt

from scitbx import linalg, matrix

from dials.algorithms.profile_model.gaussian_rs import CoordinateSystem2d
from dials.algorithms.profile_model.potato.model import (
    compute_change_of_basis_operation,
)
from dials.algorithms.profile_model.potato.util.simplex import SimpleSimplex
from dials.array_family import flex

logger = logging.getLogger("dials." + __name__)


class MarginalDistribution(object):
    """
    A class to compute useful stuff about the marginal distribution

    """

    def __init__(self, S, dS, d2S=None):

        # Compute the marginal variance
        self.S = S[8]

        # Compute the marginal derivatives
        self.dS = flex.double(d[8] for d in dS)

        # Compute the marginal second derivatives
        if d2S is not None:
            self.d2S = flex.double([d2[8] for d2 in d2S])
            self.d2S.reshape(d2S.accessor())
        else:
            self.d2S = None

    def sigma(self):
        """
        Return the marginal sigma

        """
        return self.S

    def first_derivatives(self):
        """
        Return the marginal first derivatives

        """
        return self.dS

    def second_derivatives(self):
        """
        Return the maginal second derivatives

        """
        return self.d2S


class ConditionalDistribution(object):
    """
    A class to compute useful stuff about the conditional distribution

    """

    def __init__(self, s0, s2, S, dS, d2S=None):

        self._S = S
        self._dS = dS
        self._d2S = d2S

        # Partition the covariance matrix
        S11 = matrix.sqr((S[0], S[1], S[3], S[4]))
        S12 = matrix.col((S[2], S[5]))
        S21 = matrix.col((S[6], S[7])).transpose()
        S22 = S[8]

        # The partitioned mean vector
        mu1 = matrix.col((0, 0))
        mu2 = s2.length()
        a = s0.length()

        # The epsilon
        self.epsilon = a - mu2

        # Compute the conditional mean
        self.mubar = mu1 + S12 * (1 / S22) * self.epsilon

        # Compute the conditional covariance matrix
        self.Sbar = S11 - S12 * (1 / S22) * S21

        # Set to None and compute on demand
        self.dSbar = None
        self.dmbar = None
        self.d2Sbar = None
        self.d2mbar = None

    def mean(self):
        """
        Return the conditional mean

        """
        return self.mubar

    def sigma(self):
        """
        Return the conditional sigma

        """
        return self.Sbar

    def first_derivatives_of_sigma(self):
        """
        Return the marginal first derivatives

        """

        def compute_dSbar(S, dS):

            S12 = matrix.col((S[2], S[5]))
            S21 = matrix.col((S[6], S[7])).transpose()
            S22 = S[8]

            dS11 = matrix.sqr((dS[0], dS[1], dS[3], dS[4]))
            dS12 = matrix.col((dS[2], dS[5]))
            dS21 = matrix.col((dS[6], dS[7])).transpose()
            dS22 = dS[8]

            S22_inv = 1 / S22

            A = dS11
            B = S12 * S22_inv * dS22 * S22_inv * S21
            C = S12 * S22_inv * dS21
            D = dS12 * S22_inv * S21
            return A + B - (C + D)

        if self.dSbar is None:
            self.dSbar = [compute_dSbar(self._S, d) for d in self._dS]

        return self.dSbar

    def first_derivatives_of_mean(self):
        """
        Return the marginal first derivatives

        """

        def compute_dmbar(S, dS):

            S12 = matrix.col((S[2], S[5]))
            # S21 = matrix.col((S[6], S[7])).transpose()
            S22 = S[8]

            # dS11 = matrix.sqr((dS[0], dS[1], dS[3], dS[4]))
            dS12 = matrix.col((dS[2], dS[5]))
            # dS21 = matrix.col((dS[6], dS[7])).transpose()
            dS22 = dS[8]

            S22_inv = 1 / S22

            A = dS12 * S22_inv
            B = -S12 * S22_inv * dS22 * S22_inv
            return (A + B) * self.epsilon

        if self.dmbar is None:
            self.dmbar = [compute_dmbar(self._S, d) for d in self._dS]

        return self.dmbar

    def second_derivatives_of_mean(self):
        """
        Return the maginal second derivatives

        """

        def compute_d2mbar(S, dSi, dSj, d2S):

            S12 = matrix.col((S[2], S[5]))
            # S21 = matrix.col((S[6], S[7])).transpose()
            S22 = S[8]

            dSi12 = matrix.col((dSi[2], dSi[5]))
            # dSi21 = matrix.col((dSi[6], dSi[7])).transpose()
            dSi22 = dSi[8]

            dSj12 = matrix.col((dSj[2], dSj[5]))
            # dSj21 = matrix.col((dSj[6], dSj[7])).transpose()
            dSj22 = dSj[8]

            # d2S11 = matrix.sqr((d2S[0], d2S[1], d2S[3], d2S[4]))
            d2S12 = matrix.col((d2S[2], d2S[5]))
            # d2S21 = matrix.col((d2S[6], d2S[7])).transpose()
            d2S22 = d2S[8]

            S22_inv = 1 / S22

            A = d2S12 * S22_inv
            B = dSi12 * S22_inv * dSj22 * S22_inv
            C = dSj12 * S22_inv * dSi22 * S22_inv

            D = S12 * S22_inv * dSj22 * S22_inv * dSi22 * S22_inv
            E = S12 * S22_inv * d2S22 * S22_inv
            F = S12 * S22_inv * dSi22 * S22_inv * dSj22 * S22_inv

            return (A - B - C + D - E + F) * self.epsilon

        if self.d2mbar is None:
            self.d2mbar = [
                [
                    compute_d2mbar(self._S, self._dS[i], self._dS[j], self._d2S[i, j])
                    for j in range(self._d2S.all()[1])
                ]
                for i in range(self._d2S.all()[0])
            ]

        return self.d2mbar

    def second_derivatives_of_sigma(self):
        """
        Return the maginal second derivatives

        """

        def compute_d2S(S, dSi, dSj, d2S):

            S12 = matrix.col((S[2], S[5]))
            S21 = matrix.col((S[6], S[7])).transpose()
            S22 = S[8]

            dSi12 = matrix.col((dSi[2], dSi[5]))
            dSi21 = matrix.col((dSi[6], dSi[7])).transpose()
            dSi22 = dSi[8]

            dSj12 = matrix.col((dSj[2], dSj[5]))
            dSj21 = matrix.col((dSj[6], dSj[7])).transpose()
            dSj22 = dSj[8]

            d2S11 = matrix.sqr((d2S[0], d2S[1], d2S[3], d2S[4]))
            d2S12 = matrix.col((d2S[2], d2S[5]))
            d2S21 = matrix.col((d2S[6], d2S[7])).transpose()
            d2S22 = d2S[8]

            S22_inv = 1 / S22

            A = d2S11
            B = dSj12 * S22_inv * dSi22 * S22_inv * S21
            C = S12 * S22_inv * dSj22 * S22_inv * dSi22 * S22_inv * S21

            D = S12 * S22_inv * d2S22 * S22_inv * S21
            E = S12 * S22_inv * dSi22 * S22_inv * dSj22 * S22_inv * S21
            F = S12 * S22_inv * dSi22 * S22_inv * dSj21

            G = dSj12 * S22_inv * dSi21
            H = S12 * S22_inv * dSj22 * S22_inv * dSi21
            I = S12 * S22_inv * d2S21

            J = d2S12 * S22_inv * S21
            K = dSi12 * S22_inv * (dSj22) * S22_inv * S21
            L = dSi12 * S22_inv * dSj21

            return A + B - C + D - E + F - G + H - I - J + K - L

        if self.d2Sbar is None:
            self.d2Sbar = [
                [
                    compute_d2S(self._S, self._dS[i], self._dS[j], self._d2S[i, j])
                    for j in range(self._d2S.all()[1])
                ]
                for i in range(self._d2S.all()[0])
            ]

        return self.d2Sbar


def rotate_mat3_double(R, A):
    """
    Helper function to rotate a flex.mat3_double array of matrices

    """
    accessor = A.accessor()
    RAR = flex.mat3_double([R * matrix.sqr(a) * R.transpose() for a in A])
    RAR.reshape(accessor)
    return RAR


class ReflectionData(object):
    def __init__(self, model, s0, s2, ctot, mobs, sobs, second_derivatives=False):

        # Save stuff
        self.model = model
        self.s0 = s0
        self.s2 = s2
        self.ctot = ctot
        self.mobs = mobs
        self.sobs = sobs

        # Compute the change of basis
        self.R = compute_change_of_basis_operation(s0, s2)

        # Rotate the covariance matrix
        self.S = self.R * model.sigma() * self.R.transpose()

        # Rotate the first derivative matrices
        self.dS = rotate_mat3_double(self.R, model.first_derivatives())

        # Rotate the first derivative matrices
        if second_derivatives:
            self.d2S = rotate_mat3_double(self.R, model.second_derivatives())
        else:
            self.d2S = None

        # Construct the marginal distribution
        self.marginal = MarginalDistribution(self.S, self.dS, self.d2S)

        # Construct the conditional distribution
        self.conditional = ConditionalDistribution(s0, s2, self.S, self.dS, self.d2S)

    def log_likelihood(self):
        """
        Compute the log likelihood for the reflection

        """

        # Get data
        s0 = self.s0
        s2 = self.s2
        ctot = self.ctot
        mobs = self.mobs
        Sobs = self.sobs

        # Get info about the marginal
        S22 = self.marginal.sigma()
        S22_inv = 1 / S22

        # Get info about the conditional
        Sbar = self.conditional.sigma()
        mubar = self.conditional.mean()
        Sbar_inv = Sbar.inverse()
        Sbar_det = Sbar.determinant()

        # Compute the marginal likelihood
        m_d = s0.length() - s2.length()
        m_lnL = log(S22) + S22_inv * m_d ** 2

        # Compute the conditional likelihood
        c_d = mobs - mubar
        c_lnL = ctot * (
            log(Sbar_det) + (Sbar_inv * (Sobs + c_d * c_d.transpose())).trace()
        )
        # c_lnL = ctot*(log(Sbar_det) + (Sbar_inv * Sobs).trace())

        # Return the joint likelihood
        return -0.5 * (m_lnL + c_lnL)

    def first_derivatives(self):
        """
        Compute the first derivatives

        """
        # Get data
        s0 = self.s0
        s2 = self.s2
        ctot = self.ctot
        mobs = self.mobs
        Sobs = self.sobs

        # Get info about marginal distribution
        S22 = self.marginal.sigma()
        dS22 = self.marginal.first_derivatives()
        S22_inv = 1 / S22

        # Get info about conditional distribution
        Sbar = self.conditional.sigma()
        mubar = self.conditional.mean()
        dSbar = self.conditional.first_derivatives_of_sigma()
        dmbar = self.conditional.first_derivatives_of_mean()
        Sbar_inv = Sbar.inverse()

        # The distance from the ewald sphere
        m_d = s0.length() - s2.length()
        c_d = mobs - mubar

        # Compute the derivative wrt parameter i
        dL = flex.double()
        for i in range(len(dS22)):

            I = matrix.sqr((1, 0, 0, 1))

            U = S22_inv * dS22[i] * (1 - S22_inv * m_d ** 2)
            V = (
                Sbar_inv
                * dSbar[i]
                * ctot
                * (I - Sbar_inv * (Sobs + c_d * c_d.transpose()))
            ).trace()
            W = (-2 * ctot * Sbar_inv * c_d * dmbar[i].transpose()).trace()

            dL.append(-0.5 * (U + V + W))

            # U = S22_inv*dS22[i]*(1 - S22_inv*m_d**2)
            # V = (Sbar_inv*dSbar[i]*ctot*(I - Sbar_inv*Sobs)).trace()

            # dL.append(-0.5*(U+V))

        # Return the derivative of the log likelihood
        return dL

    def second_derivatives(self):
        """
        Compute the second derivatives

        """
        # Get data
        s0 = self.s0
        s2 = self.s2
        ctot = self.ctot
        mobs = self.mobs
        Sobs = self.sobs

        # Get info about marginal distribution
        S22 = self.marginal.sigma()
        dS22 = self.marginal.first_derivatives()
        d2S22 = self.marginal.second_derivatives()
        S22_inv = 1 / S22

        # Get info about conditional distribution
        Sbar = self.conditional.sigma()
        dSbar = self.conditional.first_derivatives_of_sigma()
        d2Sbar = self.conditional.second_derivatives_of_sigma()
        mubar = self.conditional.mean()
        dmbar = self.conditional.first_derivatives_of_mean()
        d2mbar = self.conditional.second_derivatives_of_mean()
        Sbar_inv = Sbar.inverse()

        # The distance from the ewald sphere
        m_d = s0.length() - s2.length()
        c_d = mobs - mubar

        # Compute the derivative wrt parameter i j
        d2L = flex.double(d2S22.accessor())
        for j in range(d2S22.all()[0]):
            for i in range(d2S22.all()[1]):

                I = matrix.sqr((1, 0, 0, 1))

                A1 = S22_inv * d2S22[j, i] * (1 - S22_inv * m_d ** 2)
                A2 = (
                    S22_inv * dS22[j] * S22_inv * dS22[i] * (1 - 2 * S22_inv * m_d ** 2)
                )
                B1 = (
                    Sbar_inv
                    * d2Sbar[j][i]
                    * ctot
                    * (I - Sbar_inv * (Sobs + c_d * c_d.transpose()))
                )
                B2 = (
                    Sbar_inv
                    * dSbar[j]
                    * Sbar_inv
                    * dSbar[i]
                    * ctot
                    * (I - 2 * Sbar_inv * (Sobs + c_d * c_d.transpose()))
                )
                B3 = (
                    Sbar_inv
                    * dSbar[i]
                    * 2
                    * Sbar_inv
                    * ctot
                    * c_d
                    * dmbar[j].transpose()
                )
                B4 = (
                    Sbar_inv
                    * dSbar[j]
                    * 2
                    * Sbar_inv
                    * ctot
                    * c_d
                    * dmbar[i].transpose()
                )
                B5 = 2 * Sbar_inv * ctot * dmbar[j] * dmbar[i].transpose()
                B6 = 2 * Sbar_inv * ctot * c_d * d2mbar[j][i].transpose()
                U = A1 - A2
                V = (B1 - B2 + B3 + B4 + B5 - B6).trace()
                d2L[j, i] = -0.5 * (U + V)

        # Return the second derivative of the log likelihood
        return d2L

    def fisher_information(self):
        """
        Compute the fisher information

        """
        # Get data
        # s0 = self.s0
        # s2 = self.s2
        ctot = self.ctot
        # mobs = self.mobs
        # Sobs = self.sobs

        # Get info about marginal distribution
        S22 = self.marginal.sigma()
        dS22 = self.marginal.first_derivatives()
        S22_inv = 1 / S22

        # Get info about conditional distribution
        Sbar = self.conditional.sigma()
        # mbar = self.conditional.mean()
        dSbar = self.conditional.first_derivatives_of_sigma()
        dmbar = self.conditional.first_derivatives_of_mean()
        Sbar_inv = Sbar.inverse()

        # Compute the fisher information wrt parameter i j
        I = flex.double(flex.grid(len(dS22), len(dS22)))
        for j in range(len(dS22)):
            for i in range(len(dS22)):

                U = S22_inv * dS22[j] * S22_inv * dS22[i]
                V = (Sbar_inv * dSbar[j] * Sbar_inv * dSbar[i] * ctot).trace()
                W = ctot * (dmbar[i].transpose() * Sbar_inv * dmbar[j])[0]
                I[j, i] = 0.5 * (U + V) + W

                # U = S22_inv*dS22[j]*S22_inv*dS22[i]
                # V = (Sbar_inv*dSbar[j]*Sbar_inv*dSbar[i]*ctot).trace()
                # I[j,i] = 0.5*(U+V)

        return I


class MaximumLikelihoodTarget(object):
    def __init__(self, model, s0, s2_list, ctot_list, mobs_list, sobs_list):

        # Check input
        assert len(s2_list) == len(ctot_list)
        assert len(s2_list) == len(mobs_list)
        assert len(s2_list) == sobs_list.all()[0]

        # Compute the change of basis for each reflection
        self.data = []
        for i in range(len(s2_list)):
            self.data.append(
                ReflectionData(
                    model,
                    s0,
                    matrix.col(s2_list[i]),
                    ctot_list[i],
                    matrix.col(mobs_list[i]),
                    matrix.sqr(sobs_list[i : i + 1, :]),
                )
            )

    def log_likelihood(self):
        """
        The joint log likelihood

        """
        lnL = 0
        for i in range(len(self.data)):
            lnL += self.data[i].log_likelihood()
        return lnL

    def first_derivatives(self):
        """
        The joint first derivatives

        """
        dL = 0
        for i in range(len(self.data)):
            dL += self.data[i].first_derivatives()
        return dL

    def second_derivatives(self):
        """
        The joint second derivatives

        """
        d2L = 0
        for i in range(len(self.data)):
            d2L += self.data[i].second_derivatives()
        return d2L

    def fisher_information(self):
        """
        The joint fisher information

        """
        I = 0
        for i in range(len(self.data)):
            I += self.data[i].fisher_information()
        return I


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
    min_delta = min(tolerance, tolerance / p.length())
    while delta > min_delta:
        fb = func(x + delta * p)
        if fb <= fa:
            return delta
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

    def __init__(self, x0, max_iter=1000, tolerance=1e-7):
        """
        Configure the algorithm

        :param x0: The initial parameter estimates
        :param max_iter: The maximum number of iterations
        :param tolerance: The parameter tolerance

        """
        self.x0 = matrix.col(x0)
        self.max_iter = max_iter
        self.tolerance = tolerance

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
            # (e.g. 2 reflections) do an iteration of gradient descent
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

            # Update the parameter
            x0 = x

        # Save the parameters
        self.num_iter = it + 1
        self.parameters = x

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
        parameterisation,
        s0,
        s2_list,
        ctot_list,
        mobs_list,
        sobs_list,
        max_iter=1000,
        tolerance=1e-7,
    ):
        """
        Initialise the algorithm:

        """

        # Initialise the super class
        super(FisherScoringMaximumLikelihood, self).__init__(
            parameterisation.parameters(), max_iter=max_iter, tolerance=tolerance
        )

        # Save the parameterisation
        self.parameterisation = parameterisation

        # Save some stuff
        self.s0 = s0
        self.s2_list = s2_list
        self.ctot_list = ctot_list
        self.mobs_list = mobs_list
        self.sobs_list = sobs_list

        # Store the parameter history
        self.history = []

        # Print initial
        self.callback(self.parameterisation.parameters())

    def log_likelihood(self, x):
        """
        :param x: The parameter estimate
        :return: The log likelihood at x

        """
        return self.model(x).log_likelihood()

    def score(self, x):
        """
        :param x: The parameter estimate
        :return: The score at x

        """
        return self.model(x).first_derivatives()

    def score_and_fisher_information(self, x):
        """
        :param x: The parameter estimate
        :return: The score and fisher information at x

        """
        model = self.model(x)
        S = model.first_derivatives()
        I = model.fisher_information()
        return S, I

    def model(self, x):
        """
        :param x: The parameter estimate
        :return: The model

        """
        self.parameterisation.set_parameters(x)
        model = MaximumLikelihoodTarget(
            self.parameterisation,
            self.s0,
            self.s2_list,
            self.ctot_list,
            self.mobs_list,
            self.sobs_list,
        )
        return model

    def callback(self, x):
        """
        Handle and update in parameter values

        """
        self.parameterisation.set_parameters(x)
        sigma = self.parameterisation.sigma()
        lnL = self.log_likelihood(x)
        format_string = (
            "( %.3g, %.3g, %.3g, %.3g, %.3g, %.3g, %.3g, %.3g, %.3g ): L = %f"
        )
        logger.info(format_string % (tuple(sigma) + (lnL,)))
        self.history.append(x)


class ProfileRefiner(object):
    """
    High level profile refiner class that handles book keeping etc

    """

    def __init__(self, parameterisation, data):
        """
        Set the data and initial parameters

        """
        self.s0 = data.s0
        self.s2_list = data.s2_list
        self.ctot_list = data.ctot_list
        self.mobs_list = data.mobs_list
        self.sobs_list = data.sobs_list
        self.parameterisation = parameterisation

    def refine(self):
        """
        Perform the profile refinement

        """
        if False:
            self.refine_simplex()
        else:
            self.refine_fisher_scoring()

    def refine_simplex(self):
        """
        Perform the profile refinement

        """

        class Target(object):
            def __init__(self, model, s0, s2_list, ctot_list, xbar_list, sobs_list):
                self.model = model
                self.s0 = s0
                self.s2_list = s2_list
                self.ctot_list = ctot_list
                self.xbar_list = xbar_list
                self.sobs_list = sobs_list

            def target(self, params):
                self.model.set_parameters(params)
                ml = MaximumLikelihoodTarget(
                    self.model,
                    self.s0,
                    self.s2_list,
                    self.ctot_list,
                    self.xbar_list,
                    self.sobs_list,
                )
                lnL = ml.log_likelihood()
                sigma = self.model.sigma()
                format_string = (
                    "( %.3g, %.3g, %.3g, %.3g, %.3g, %.3g, %.3g, %.3g, %.3g ): L = %f"
                )
                logger.info(format_string % (tuple(sigma) + (lnL,)))
                return -lnL

        # Starting values for simplex
        values = flex.double(self.parameterisation.parameters())
        offset = flex.double([sqrt(1e-7) for v in values])

        # Do the simplex optimization
        optimizer = SimpleSimplex(
            values,
            offset,
            Target(
                self.parameterisation,
                self.s0,
                self.s2_list,
                self.ctot_list,
                self.mobs_list,
                self.sobs_list,
            ),
            2000,
        )

        # Get the parameters
        self.parameters = optimizer.get_solution()

        # set the parameters
        self.parameterisation.set_parameters(self.parameters)

        # Print the eigen values and vectors
        print_eigen_values_and_vectors(self.parameterisation.sigma())

    def refine_fisher_scoring(self):
        """
        Perform the profile refinement

        """

        # Initialise the algorithm
        ml = FisherScoringMaximumLikelihood(
            self.parameterisation,
            self.s0,
            self.s2_list,
            self.ctot_list,
            self.mobs_list,
            self.sobs_list,
        )

        # Solve the maximum likelihood equations
        ml.solve()

        # Get the parameters
        self.parameters = ml.parameters

        # set the parameters
        self.parameterisation.set_parameters(self.parameters)

        # Print the eigen values and vectors
        print_eigen_values_and_vectors(self.parameterisation.sigma())

        # Return the optimizer
        return ml


class ProfileRefinerData(object):
    """
    A class for holding the data needed for the profile refinement

    """

    def __init__(self, s0, s2_list, ctot_list, mobs_list, sobs_list):
        """
        Init the data

        """
        self.s0 = s0
        self.s2_list = s2_list
        self.ctot_list = ctot_list
        self.mobs_list = mobs_list
        self.sobs_list = sobs_list

    @classmethod
    def from_reflections(self, experiment, reflections):
        """
        Generate the required data from the reflections

        """

        # Get the beam vector
        s0 = matrix.col(experiment.beam.get_s0())

        # Get the reciprocal lattice vector
        s2_list = reflections["s2"]

        # Initialise the list of observed intensities and covariances
        ctot_list = flex.double(len(s2_list))
        mobs_list = flex.vec2_double(len(s2_list))
        Sobs_list = flex.double(flex.grid(len(s2_list), 4))
        Bmean = matrix.sqr((0, 0, 0, 0))

        logger.info(
            "Computing observed covariance for %d reflections" % len(reflections)
        )
        s0_length = s0.length()
        assert len(experiment.detector) == 1
        panel = experiment.detector[0]
        sbox = reflections["shoebox"]
        for r in range(len(reflections)):

            # Create the coordinate system
            cs = CoordinateSystem2d(s0, s2_list[r])

            # Get data and compute total counts
            data = sbox[r].data
            mask = sbox[r].mask
            bgrd = sbox[r].background

            # Get array of vectors
            i0 = sbox[r].bbox[0]
            j0 = sbox[r].bbox[2]
            assert data.all()[0] == 1
            X = flex.vec2_double(flex.grid(data.all()[1], data.all()[2]))
            ctot = 0
            C = flex.double(X.accessor())

            for j in range(data.all()[1]):
                for i in range(data.all()[2]):
                    c = data[0, j, i] - bgrd[0, j, i]
                    if mask[0, j, i] == 5 and c > 0:
                        ctot += c
                        ii = i + i0
                        jj = j + j0
                        s = panel.get_pixel_lab_coord((ii + 0.5, jj + 0.5))
                        s = matrix.col(s).normalize() * s0_length
                        X[j, i] = cs.from_beam_vector(s)
                        C[j, i] = c

            # Check we have a sensible number of counts
            assert ctot > 0, "BUG: strong spots should have more than 0 counts!"

            # Compute the mean vector
            xbar = matrix.col((0, 0))
            for j in range(X.all()[0]):
                for i in range(X.all()[1]):
                    x = matrix.col(X[j, i])
                    xbar += C[j, i] * x
            xbar /= ctot

            # Compute the covariance matrix
            Sobs = matrix.sqr((0, 0, 0, 0))
            for j in range(X.all()[0]):
                for i in range(X.all()[1]):
                    x = matrix.col(X[j, i])
                    Sobs += (x - xbar) * (x - xbar).transpose() * C[j, i]
            Sobs /= ctot

            # Compute the bias
            zero = matrix.col((0, 0))
            Bias_sq = (xbar - zero) * (xbar - zero).transpose()
            Bmean += Bias_sq

            # Add to the lists
            ctot_list[r] = ctot
            mobs_list[r] = xbar
            Sobs_list[r, 0] = Sobs[0]
            Sobs_list[r, 1] = Sobs[1]
            Sobs_list[r, 2] = Sobs[2]
            Sobs_list[r, 3] = Sobs[3]

        # Print some information
        logger.info(
            "I_min = %.2f, I_max = %.2f" % (flex.min(ctot_list), flex.max(ctot_list))
        )

        # Print the mean covariance
        Smean = matrix.sqr((0, 0, 0, 0))
        for r in range(Sobs_list.all()[0]):
            Smean += matrix.sqr(tuple(Sobs_list[r : r + 1, :]))
        Smean /= Sobs_list.all()[0]
        Bmean /= len(reflections)

        logger.info("")
        logger.info("Mean observed covariance:")
        print_matrix(Smean)
        print_eigen_values_and_vectors_of_observed_covariance(Smean)
        logger.info("")
        logger.info("Mean observed bias^2:")
        print_matrix(Bmean)

        # Compute the distance from the Ewald sphere
        epsilon = flex.double(
            s0.length() - matrix.col(s).length() for s in reflections["s2"]
        )
        mv = flex.mean_and_variance(epsilon)
        logger.info("")
        logger.info("Mean distance from Ewald sphere: %.3g" % mv.mean())
        logger.info(
            "Variance in distance from Ewald sphere: %.3g"
            % mv.unweighted_sample_variance()
        )

        # Return the profile refiner data
        return ProfileRefinerData(s0, s2_list, ctot_list, mobs_list, Sobs_list)


def print_eigen_values_and_vectors_of_observed_covariance(A):
    """
    Print the eigen values and vectors of a matrix

    """

    # Compute the eigen decomposition of the covariance matrix
    eigen_decomposition = linalg.eigensystem.real_symmetric(A.as_flex_double_matrix())
    Q = matrix.sqr(eigen_decomposition.vectors())
    L = matrix.diag(eigen_decomposition.values())

    # Print the matrix eigen values
    logger.info("")
    logger.info("Eigen Values:")
    logger.info("")
    print_matrix(L, indent=2)
    logger.info("")

    logger.info("Eigen Vectors:")
    logger.info("")
    print_matrix(Q, indent=2)
    logger.info("")

    logger.info("Observed covariance in degrees equivalent units")
    logger.info("C1: %.5f degrees" % (sqrt(L[0]) * 180.0 / pi))
    logger.info("C2: %.5f degrees" % (sqrt(L[3]) * 180.0 / pi))


def print_eigen_values_and_vectors(A):
    """
    Print the eigen values and vectors of a matrix

    """

    # Compute the eigen decomposition of the covariance matrix
    eigen_decomposition = linalg.eigensystem.real_symmetric(A.as_flex_double_matrix())
    Q = matrix.sqr(eigen_decomposition.vectors())
    L = matrix.diag(eigen_decomposition.values())

    # Print the matrix eigen values
    logger.info("")
    logger.info("Eigen Values:")
    logger.info("")
    print_matrix(L, indent=2)
    logger.info("")

    logger.info("Eigen Vectors:")
    logger.info("")
    print_matrix(Q, indent=2)
    logger.info("")

    logger.info("Mosaicity in degrees equivalent units")
    logger.info("M1: %.5f degrees" % (sqrt(L[0]) * 180.0 / pi))
    logger.info("M2: %.5f degrees" % (sqrt(L[4]) * 180.0 / pi))
    logger.info("M3: %.5f degrees" % (sqrt(L[8]) * 180.0 / pi))


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
    logger.info("\n".join(lines))
