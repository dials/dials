from __future__ import division

import logging
import textwrap
from math import log, pi, sqrt

from scitbx import linalg, matrix

from dials.algorithms.profile_model.potato.model import (
    compute_change_of_basis_operation,
)
from dials.algorithms.profile_model.potato.parameterisation import ReflectionModelState
from dials.algorithms.profile_model.potato.util.simplex import SimpleSimplex
from dials.array_family import flex
from dials.util import tabulate

logger = logging.getLogger("dials." + __name__)


class ConditionalDistribution(object):
    """
    A class to compute useful stuff about the conditional distribution

    """

    def __init__(self, s0, mu, dmu, S, dS):

        self._mu = mu
        self._dmu = dmu
        self._S = S
        self._dS = dS

        # Partition the covariance matrix
        S11 = matrix.sqr((S[0], S[1], S[3], S[4]))
        S12 = matrix.col((S[2], S[5]))
        S21 = matrix.col((S[6], S[7])).transpose()
        S22 = S[8]

        # The partitioned mean vector
        mu1 = matrix.col((mu[0], mu[1]))
        mu2 = mu[2]
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

        def compute_dmbar(i):

            S = self._S
            dS = self._dS[i]
            # mu = self._mu
            dmu = self._dmu[i]

            S12 = matrix.col((S[2], S[5]))
            # S21 = matrix.col((S[6], S[7])).transpose()
            S22 = S[8]

            # dS11 = matrix.sqr((dS[0], dS[1], dS[3], dS[4]))
            dS12 = matrix.col((dS[2], dS[5]))
            # dS21 = matrix.col((dS[6], dS[7])).transpose()
            dS22 = dS[8]

            S22_inv = 1 / S22

            dmu1 = matrix.col((dmu[0], dmu[1]))
            dmu2 = dmu[2]
            dep = -dmu2

            A = dmu1
            B = dS12 * S22_inv * self.epsilon
            C = -S12 * S22_inv * dS22 * S22_inv * self.epsilon
            D = S12 * S22_inv * dep
            return A + B + C + D

        if self.dmbar is None:
            self.dmbar = [compute_dmbar(i) for i in range(len(self._dS))]

        return self.dmbar


def rotate_vec3_double(R, A):
    """
    Helper function to rotate a flex.mat3_double array of matrices

    """
    accessor = A.accessor()
    RA = flex.vec3_double([R * matrix.col(a) for a in A])
    RA.reshape(accessor)
    return RA


def rotate_mat3_double(R, A):
    """
    Helper function to rotate a flex.mat3_double array of matrices

    """
    accessor = A.accessor()
    RAR = flex.mat3_double([R * matrix.sqr(a) * R.transpose() for a in A])
    RAR.reshape(accessor)
    return RAR


class ReflectionLikelihood(object):
    def __init__(self, model, s0, sp, h, ctot, mobs, sobs):

        # Save stuff
        self.model = ReflectionModelState(model, s0, h)
        self.s0 = s0
        self.sp = sp
        self.h = h
        self.ctot = ctot
        self.mobs = mobs
        self.sobs = sobs

        # Compute the change of basis
        self.R = compute_change_of_basis_operation(s0, sp)

        # The s2 vector
        self.r = self.model.get_r()
        self.s2 = s0 + self.r

        # Rotate the mean vector
        self.mu = self.R * self.s2

        # Rotate the covariance matrix
        self.S = self.R * self.model.get_sigma() * self.R.transpose()

        # Rotate the first derivative matrices
        self.dS = rotate_mat3_double(self.R, self.model.get_dS_dp())

        # Rotate the first derivative of s2
        self.dmu = rotate_vec3_double(self.R, self.model.get_dr_dp())

        # Construct the conditional distribution
        self.conditional = ConditionalDistribution(
            s0, self.mu, self.dmu, self.S, self.dS
        )

    def log_likelihood(self):
        """
        Compute the log likelihood for the reflection

        """

        # Get data
        s0 = self.s0
        # s2 = self.s2
        ctot = self.ctot
        mobs = self.mobs
        Sobs = self.sobs

        # Get info about the marginal
        S22 = self.S[8]
        S22_inv = 1 / S22
        mu2 = self.mu[2]

        # Get info about the conditional
        Sbar = self.conditional.sigma()
        mubar = self.conditional.mean()
        Sbar_inv = Sbar.inverse()
        Sbar_det = Sbar.determinant()

        # Weights for marginal and conditional components
        m_w = ctot
        c_w = ctot

        # Compute the marginal likelihood
        m_d = s0.length() - mu2
        m_lnL = m_w * (log(S22) + S22_inv * m_d ** 2)

        # Compute the conditional likelihood
        c_d = mobs - mubar
        c_lnL = c_w * (
            log(Sbar_det) + (Sbar_inv * (Sobs + c_d * c_d.transpose())).trace()
        )

        # Return the joint likelihood
        return -0.5 * (m_lnL + c_lnL)

    def first_derivatives(self):
        """
        Compute the first derivatives

        """
        # Get data
        s0 = self.s0
        # s2 = self.s2
        ctot = self.ctot
        mobs = self.mobs
        Sobs = self.sobs

        # Get info about marginal distribution
        S22 = self.S[8]
        dS22 = flex.double(self.dS[i][8] for i in range(len(self.dS)))
        S22_inv = 1 / S22
        mu2 = self.mu[2]

        # Get info about conditional distribution
        Sbar = self.conditional.sigma()
        mubar = self.conditional.mean()
        dSbar = self.conditional.first_derivatives_of_sigma()
        dmbar = self.conditional.first_derivatives_of_mean()
        Sbar_inv = Sbar.inverse()

        # The distance from the ewald sphere
        epsilon = s0.length() - mu2
        c_d = mobs - mubar

        # Weights for marginal and conditional components
        m_w = ctot
        c_w = ctot

        # Compute the derivative wrt parameter i
        dL = flex.double()
        for i in range(len(dS22)):

            dmu = self.dmu[i]
            # dmu1 = matrix.col((dmu[0], dmu[1]))
            dmu2 = dmu[2]
            dep = -dmu2

            I = matrix.sqr((1, 0, 0, 1))

            U = m_w * (
                S22_inv * dS22[i] * (1 - S22_inv * epsilon ** 2)
                + 2 * S22_inv * epsilon * dep
            )
            V = (
                c_w
                * (
                    Sbar_inv
                    * dSbar[i]
                    * (I - Sbar_inv * (Sobs + c_d * c_d.transpose()))
                ).trace()
            )
            W = c_w * (-2 * Sbar_inv * c_d * dmbar[i].transpose()).trace()
            dL.append(-0.5 * (U + V + W))

        # Return the derivative of the log likelihood
        return dL

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
        S22 = self.S[8]
        dS22 = flex.double(self.dS[i][8] for i in range(len(self.dS)))
        S22_inv = 1 / S22

        # Get info about conditional distribution
        Sbar = self.conditional.sigma()
        # mbar = self.conditional.mean()
        dSbar = self.conditional.first_derivatives_of_sigma()
        dmbar = self.conditional.first_derivatives_of_mean()
        Sbar_inv = Sbar.inverse()
        dmu = self.dmu

        # Weights for marginal and conditional components
        m_w = ctot
        c_w = ctot

        # Compute the fisher information wrt parameter i j
        I = flex.double(flex.grid(len(dS22), len(dS22)))
        for j in range(len(dS22)):
            for i in range(len(dS22)):
                U = S22_inv * dS22[j] * S22_inv * dS22[i]
                V = (Sbar_inv * dSbar[j] * Sbar_inv * dSbar[i]).trace()
                W = 2 * (Sbar_inv * dmbar[i] * dmbar[j].transpose()).trace()
                X = 2 * dmu[i][2] * S22_inv * dmu[j][2]
                I[j, i] = 0.5 * c_w * (V + W) + 0.5 * m_w * (U + X)

        return I


class MaximumLikelihoodTarget(object):
    def __init__(self, model, s0, sp_list, h_list, ctot_list, mobs_list, sobs_list):

        # Check input
        assert len(h_list) == len(sp_list)
        assert len(h_list) == len(ctot_list)
        assert len(h_list) == len(mobs_list)
        assert len(h_list) == sobs_list.all()[0]

        # Save the model
        self.model = model

        # Compute the change of basis for each reflection
        self.data = []
        for i in range(len(h_list)):
            self.data.append(
                ReflectionLikelihood(
                    model,
                    s0,
                    matrix.col(sp_list[i]),
                    matrix.col(h_list[i]),
                    ctot_list[i],
                    matrix.col(mobs_list[i]),
                    matrix.sqr(sobs_list[i : i + 1, :]),
                )
            )

    def mse(self):
        """
        The MSE in local reflection coordinates

        """
        mse = 0
        for i in range(len(self.data)):
            mbar = self.data[i].conditional.mean()
            xobs = self.data[i].mobs
            mse += (xobs - mbar).dot(xobs - mbar)
        mse /= len(self.data)
        return mse

    def rmsd(self):
        """
        The RMSD in pixels

        """
        mse = matrix.col((0, 0))
        for i in range(len(self.data)):
            R = self.data[i].R
            mbar = self.data[i].conditional.mean()
            xobs = self.data[i].mobs
            s0 = self.data[i].s0
            s1 = R.transpose() * matrix.col((mbar[0], mbar[1], s0.length()))
            s3 = R.transpose() * matrix.col((xobs[0], xobs[1], s0.length()))
            xyzcal = self.model.experiment.detector[0].get_ray_intersection_px(s1)
            xyzobs = self.model.experiment.detector[0].get_ray_intersection_px(s3)
            r_x = xyzcal[0] - xyzobs[0]
            r_y = xyzcal[1] - xyzobs[1]
            mse += matrix.col((r_x ** 2, r_y ** 2))
        mse /= len(self.data)
        return matrix.col((sqrt(mse[0]), sqrt(mse[1])))

    def log_likelihood(self):
        """
        The joint log likelihood

        """
        lnL = 0
        for i in range(len(self.data)):
            lnL += self.data[i].log_likelihood()
        return lnL

    def jacobian(self):
        """
        Return the Jacobean

        """
        J = []
        for i in range(len(self.data)):
            J.append(list(self.data[i].first_derivatives()))
        return flex.double(J)

    def first_derivatives(self):
        """
        The joint first derivatives

        """
        dL = 0
        for i in range(len(self.data)):
            dL += self.data[i].first_derivatives()
        return dL

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
        model,
        s0,
        sp_list,
        h_list,
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
            model.get_active_parameters(), max_iter=max_iter, tolerance=tolerance
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

        # Store the parameter history
        self.history = []

        # Print initial
        self.callback(self.model.get_active_parameters())

    def log_likelihood(self, x):
        """
        :param x: The parameter estimate
        :return: The log likelihood at x

        """
        return self.target(x).log_likelihood()

    def score(self, x):
        """
        :param x: The parameter estimate
        :return: The score at x

        """
        return self.target(x).first_derivatives()

    def score_and_fisher_information(self, x):
        """
        :param x: The parameter estimate
        :return: The score and fisher information at x

        """
        model = self.target(x)
        S = model.first_derivatives()
        I = model.fisher_information()
        return S, I

    def mse(self, x):
        """
        :param x: The parameter estimate
        :return: The MSE at x

        """
        return self.target(x).mse()

    def rmsd(self, x):
        """
        :param x: The parameter estimate
        :return: The RMSD at x

        """
        return self.target(x).rmsd()

    def jacobian(self, x):
        """
        :param x: The parameter estimate
        :return: The Jacobian at x

        """
        return self.target(x).jacobian()

    def target(self, x):
        """
        :param x: The parameter estimate
        :return: The model

        """
        self.model.set_active_parameters(x)
        target = MaximumLikelihoodTarget(
            self.model,
            self.s0,
            self.sp_list,
            self.h_list,
            self.ctot_list,
            self.mobs_list,
            self.sobs_list,
        )
        return target

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
        self.model.set_active_parameters(x)
        lnL = self.log_likelihood(x)
        mse = self.mse(x)
        rmsd = self.rmsd(x)

        # Get the unit cell
        unit_cell = self.model.get_unit_cell().parameters()

        # Get some matrices
        U = self.model.get_U()
        M = self.model.get_M()

        # Print some information
        format_string1 = "  Unit cell: (%.3f, %.3f, %.3f, %.3f, %.3f, %.3f)"
        format_string2 = "  | % .6f % .6f % .6f |"
        format_string3 = "  | % .2e % .2e % .2e |"

        logger.info(f"\nIteration: {len(self.history)}")
        if not self.model.is_unit_cell_fixed():
            logger.info("\n" + format_string1 % unit_cell)
        if not self.model.is_orientation_fixed():
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
        if not self.model.is_mosaic_spread_fixed():
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

    def __init__(self, state, data):
        """
        Set the data and initial parameters

        """
        self.s0 = data.s0
        self.h_list = data.h_list
        self.sp_list = data.sp_list
        self.ctot_list = data.ctot_list
        self.mobs_list = data.mobs_list
        self.sobs_list = data.sobs_list
        self.state = state
        self.history = []

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
            def __init__(
                self, state, s0, sp_list, h_list, ctot_list, xbar_list, sobs_list
            ):
                self.state = state
                self.s0 = s0
                self.sp_list = sp_list
                self.h_list = h_list
                self.ctot_list = ctot_list
                self.xbar_list = xbar_list
                self.sobs_list = sobs_list

            def target(self, params):
                self.state.set_active_parameters(params)
                ml = MaximumLikelihoodTarget(
                    self.state,
                    self.s0,
                    self.sp_list,
                    self.h_list,
                    self.ctot_list,
                    self.xbar_list,
                    self.sobs_list,
                )
                lnL = ml.log_likelihood()
                sigma = self.state.get_M()
                format_string = (
                    "( %.3g, %.3g, %.3g, %.3g, %.3g, %.3g, %.3g, %.3g, %.3g ): L = %f"
                )
                logger.info(format_string % (tuple(sigma) + (lnL,)))
                return -lnL

        # Starting values for simplex
        values = flex.double(self.state.get_active_parameters())
        offset = flex.double([sqrt(1e-7) for v in values])

        # Do the simplex optimization
        optimizer = SimpleSimplex(
            values,
            offset,
            Target(
                self.state,
                self.s0,
                self.sp_list,
                self.h_list,
                self.ctot_list,
                self.mobs_list,
                self.sobs_list,
            ),
            2000,
        )

        # Get the parameters
        self.parameters = optimizer.get_solution()

        # set the parameters
        self.state.set_active_parameters(self.parameters)

        # Print the eigen values and vectors
        print_eigen_values_and_vectors(self.state.get_M())

    def refine_fisher_scoring(self):
        """
        Perform the profile refinement

        """

        # Print information
        logger.info("\nComponents to refine:")
        logger.info(" Orientation:       %s" % (not self.state.is_orientation_fixed()))
        logger.info(" Unit cell:         %s" % (not self.state.is_unit_cell_fixed()))
        logger.info(
            " RLP mosaicity:     %s" % (not self.state.is_mosaic_spread_fixed())
        )
        logger.info(
            " Wavelength spread: %s\n" % (not self.state.is_wavelength_spread_fixed())
        )

        # Initialise the algorithm
        ml = FisherScoringMaximumLikelihood(
            self.state,
            self.s0,
            self.sp_list,
            self.h_list,
            self.ctot_list,
            self.mobs_list,
            self.sobs_list,
        )

        # Solve the maximum likelihood equations
        ml.solve()

        # Get the parameters
        self.parameters = flex.double(ml.parameters)

        # set the parameters
        self.state.set_active_parameters(self.parameters)

        # Print summary table of refinement.
        rows = []
        headers = ["Iteration", "likelihood", "RMSD (pixel) X,Y"]
        for i, h in enumerate(ml.history):
            l = h["likelihood"]
            rmsd = h["rmsd"]
            rows.append([str(i), f"{l:.4f}", f"{rmsd[0]:.3f}, {rmsd[1]:.3f}"])
        logger.info(
            "\nRefinement steps:\n\n" + textwrap.indent(tabulate(rows, headers), " ")
        )

        # Print the eigen values and vectors of sigma_m
        if not self.state.is_mosaic_spread_fixed():
            logger.info("\nDecomposition of Sigma_M:")
            print_eigen_values_and_vectors(self.state.get_M())

        # Save the history
        self.history = ml.history

        # Return the optimizer
        return ml

    def correlation(self):
        """
        Return the correlation matrix between parameters

        """
        # Initialise the algorithm
        ml = FisherScoringMaximumLikelihood(
            self.state,
            self.s0,
            self.sp_list,
            self.h_list,
            self.ctot_list,
            self.mobs_list,
            self.sobs_list,
        )
        return ml.correlation(self.state.get_active_parameters())

    def labels(self):
        """
        Return parameter labels

        """
        return self.state.get_labels()


class RefinerData(object):
    """
    A class for holding the data needed for the profile refinement

    """

    def __init__(self, s0, sp_list, h_list, ctot_list, mobs_list, sobs_list):
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

    @classmethod
    def from_reflections(self, experiment, reflections):
        """
        Generate the required data from the reflections

        """

        # Get the beam vector
        s0 = matrix.col(experiment.beam.get_s0())

        # Get the reciprocal lattice vector
        h_list = reflections["miller_index"]

        # Initialise the list of observed intensities and covariances
        sp_list = flex.vec3_double(len(h_list))
        ctot_list = flex.double(len(h_list))
        mobs_list = flex.vec2_double(len(h_list))
        Sobs_list = flex.double(flex.grid(len(h_list), 4))
        # Bmean = matrix.sqr((0, 0, 0, 0))

        # SSS = 0
        logger.info(
            "Computing observed covariance for %d reflections" % len(reflections)
        )
        s0_length = s0.length()
        assert len(experiment.detector) == 1
        panel = experiment.detector[0]
        sbox = reflections["shoebox"]
        xyzobs = reflections["xyzobs.px.value"]
        for r in range(len(reflections)):

            # The vector to the pixel centroid
            sp = (
                matrix.col(panel.get_pixel_lab_coord(xyzobs[r][0:2])).normalize()
                * s0.length()
            )

            # Compute change of basis
            R = compute_change_of_basis_operation(s0, sp)

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
                    # if mask[0,j,i] & (1 | 4 | 8) == (1 | 4 | 8) and c > 0:
                    if mask[0, j, i] & (1 | 4) == (1 | 4) and c > 0:
                        ctot += c
                        ii = i + i0
                        jj = j + j0
                        s = panel.get_pixel_lab_coord((ii + 0.5, jj + 0.5))
                        s = matrix.col(s).normalize() * s0_length
                        e = R * s
                        X[j, i] = (e[0], e[1])
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
            assert Sobs[0] > 0, "BUG: variance must be > 0"
            assert Sobs[3] > 0, "BUG: variance must be > 0"

            # Compute the bias
            # zero = matrix.col((0, 0))
            # Bias_sq = (xbar - zero)*(xbar - zero).transpose()
            # Bmean += Bias_sq

            # ctot += 10000

            # Add to the lists
            sp_list[r] = sp
            ctot_list[r] = ctot
            mobs_list[r] = xbar
            Sobs_list[r, 0] = Sobs[0]
            Sobs_list[r, 1] = Sobs[1]
            Sobs_list[r, 2] = Sobs[2]
            Sobs_list[r, 3] = Sobs[3]

        # Print some information
        logger.info("")
        logger.info(
            "I_min = %.2f, I_max = %.2f" % (flex.min(ctot_list), flex.max(ctot_list))
        )

        # Sometimes a single reflection might have an enormouse intensity for
        # whatever reason and since we weight by intensity, this can cause the
        # refinement to be dominated by these reflections. Therefore, if the
        # intensity is greater than some value, damp the weighting accordingly
        def damp_outlier_intensity_weights(ctot_list):
            n = ctot_list.size()
            sorted_ctot = sorted(ctot_list)
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
        Smean = matrix.sqr((0, 0, 0, 0))
        for r in range(Sobs_list.all()[0]):
            Smean += matrix.sqr(tuple(Sobs_list[r : r + 1, :]))
        Smean /= Sobs_list.all()[0]
        # Bmean /= len(reflections)

        logger.info("")
        logger.info("Mean observed covariance:")
        print_matrix(Smean)
        print_eigen_values_and_vectors_of_observed_covariance(Smean, s0)
        # logger.info("")
        # logger.info("Mean observed bias^2:")
        # print_matrix(Bmean)

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
        return RefinerData(s0, sp_list, h_list, ctot_list, mobs_list, Sobs_list)


def print_eigen_values_and_vectors_of_observed_covariance(A, s0):
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
    logger.info("C1: %.5f degrees" % (sqrt(L[0]) * (180.0 / pi) / s0.length()))
    logger.info("C2: %.5f degrees" % (sqrt(L[3]) * (180.0 / pi) / s0.length()))


def print_eigen_values_and_vectors(A):
    """
    Print the eigen values and vectors of a matrix

    """

    # Compute the eigen decomposition of the covariance matrix
    eigen_decomposition = linalg.eigensystem.real_symmetric(A.as_flex_double_matrix())
    Q = matrix.sqr(eigen_decomposition.vectors())
    L = matrix.diag(eigen_decomposition.values())

    # Print the matrix eigen values
    logger.info(" ")
    logger.info(" Eigen Values:")
    logger.info(" ")
    print_matrix(L, indent=2)
    logger.info(" ")

    logger.info(" Eigen Vectors:")
    logger.info(" ")
    print_matrix(Q, indent=2)
    logger.info(" ")

    logger.info(" Mosaicity in degrees equivalent units")
    logger.info(" M1: %.5f degrees" % (sqrt(L[0]) * 180.0 / pi))
    logger.info(" M2: %.5f degrees" % (sqrt(L[4]) * 180.0 / pi))
    logger.info(" M3: %.5f degrees" % (sqrt(L[8]) * 180.0 / pi))


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
