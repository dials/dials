from __future__ import division

from typing import Optional

import numpy as np

from scitbx import matrix

from dials.algorithms.refinement.parameterisation.crystal_parameters import (
    CrystalOrientationParameterisation,
    CrystalUnitCellParameterisation,
)
from dials.array_family import flex


class BaseParameterisation(object):
    def __init__(self, params: Optional[np.array] = None) -> None:
        """
        Initialise with the parameters

        """
        if params is not None:
            assert len(params) == self.num_parameters()
            self.params = params
        else:
            self.params = np.array([0.0] * self.num_parameters(), dtype=np.float)

    @property
    def parameters(self):
        """
        Return the parameters

        """
        return self.params

    @parameters.setter
    def parameters(self, params: np.array) -> None:
        assert len(params) == self.num_parameters()
        self.params = params


class Simple1MosaicityParameterisation(BaseParameterisation):
    """
    A simple mosaicity parameterisation that uses 1 parameter to describe a
    multivariate normal reciprocal lattice profile. Sigma is enforced as positive
    definite by parameterising using the cholesky decomposition.

    M = | b1 0  0  |
        |  0 b1 0  |
        |  0  0 b1 |

    S = M*M^T

    """

    @staticmethod
    def is_angular():
        """
        Is angular

        """
        return False

    @staticmethod
    def num_parameters():
        """
        Get the number of parameters

        """
        return 1

    def sigma(self):
        """
        Compute the covariance matrix of the MVN from the parameters

        """
        M = matrix.sqr(
            (self.params[0], 0, 0, 0, self.params[0], 0, 0, 0, self.params[0])
        )
        return M * M.transpose()

    def first_derivatives(self):
        """
        Compute the first derivatives of Sigma w.r.t the parameters

        """
        b1 = self.params[0]

        dSdb1 = (2 * b1, 0, 0, 0, 2 * b1, 0, 0, 0, 2 * b1)

        return flex.mat3_double([dSdb1])


class Simple6MosaicityParameterisation(BaseParameterisation):
    """
    A simple mosaicity parameterisation that uses 6 parameters to describe a
    multivariate normal reciprocal lattice profile. Sigma is enforced as positive
    definite by parameterising using the cholesky decomposition.

    M = | b1 0  0  |
        | b2 b3 0  |
        | b4 b5 b6 |

    S = M*M^T

    """

    @staticmethod
    def is_angular():
        """
        Is angular

        """
        return False

    @staticmethod
    def num_parameters():
        """
        Get the number of parameters

        """
        return 6

    def sigma(self):
        """
        Compute the covariance matrix of the MVN from the parameters

        """
        M = matrix.sqr(
            (
                self.params[0],
                0,
                0,
                self.params[1],
                self.params[2],
                0,
                self.params[3],
                self.params[4],
                self.params[5],
            )
        )
        return M * M.transpose()

    def first_derivatives(self):
        """
        Compute the first derivatives of Sigma w.r.t the parameters

        """
        b1, b2, b3, b4, b5, b6 = self.params

        dSdb1 = (2 * b1, b2, b4, b2, 0, 0, b4, 0, 0)

        dSdb2 = (0, b1, 0, b1, 2 * b2, b4, 0, b4, 0)

        dSdb3 = (0, 0, 0, 0, 2 * b3, b5, 0, b5, 0)

        dSdb4 = (0, 0, b1, 0, 0, b2, b1, b2, 2 * b4)

        dSdb5 = (0, 0, 0, 0, 0, b3, 0, b3, 2 * b5)

        dSdb6 = (0, 0, 0, 0, 0, 0, 0, 0, 2 * b6)

        return flex.mat3_double([dSdb1, dSdb2, dSdb3, dSdb4, dSdb5, dSdb6])

    def second_derivatives(self):
        """
        Compute the second derivatives of Sigma w.r.t the parameters

        """
        b1, b2, b3, b4, b5, b6 = self.params

        zero = (0, 0, 0, 0, 0, 0, 0, 0, 0)

        d11 = (2, 0, 0, 0, 0, 0, 0, 0, 0)
        d12 = (0, 1, 0, 1, 0, 0, 0, 0, 0)
        d13 = zero
        d14 = (0, 0, 1, 0, 0, 0, 1, 0, 0)
        d15 = zero
        d16 = zero

        d21 = d12
        d22 = (0, 0, 0, 0, 2, 0, 0, 0, 0)
        d23 = zero
        d24 = (0, 0, 0, 0, 0, 1, 0, 1, 0)
        d25 = zero
        d26 = zero

        d31 = zero
        d32 = zero
        d33 = (0, 0, 0, 0, 2, 0, 0, 0, 0)
        d34 = zero
        d35 = (0, 0, 0, 0, 0, 1, 0, 1, 0)
        d36 = zero

        d41 = d14
        d42 = d24
        d43 = zero
        d44 = (0, 0, 0, 0, 0, 0, 0, 0, 2)
        d45 = zero
        d46 = zero

        d51 = zero
        d52 = zero
        d53 = d35
        d54 = zero
        d55 = (0, 0, 0, 0, 0, 0, 0, 0, 2)
        d56 = zero

        d61 = zero
        d62 = zero
        d63 = zero
        d64 = zero
        d65 = zero
        d66 = (0, 0, 0, 0, 0, 0, 0, 0, 2)

        d2 = flex.mat3_double(
            [
                d11,
                d12,
                d13,
                d14,
                d15,
                d16,
                d21,
                d22,
                d23,
                d24,
                d25,
                d26,
                d31,
                d32,
                d33,
                d34,
                d35,
                d36,
                d41,
                d42,
                d43,
                d44,
                d45,
                d46,
                d51,
                d52,
                d53,
                d54,
                d55,
                d56,
                d61,
                d62,
                d63,
                d64,
                d65,
                d66,
            ]
        )
        d2.reshape(flex.grid(6, 6))

        return d2


class WavelengthSpreadParameterisation(BaseParameterisation):
    """
    A simple mosaicity parameterisation that uses 1 parameters to describe a
    multivariate normal wavelength spread. Sigma is enforced as positive
    definite by parameterising using the cholesky decomposition.

    L = | 0 0 0  |
        | 0 0 0  |
        | 0 0 l1 |

    S = L*L^T

    """

    @staticmethod
    def num_parameters():
        """
        Get the number of parameters

        """
        return 1

    def sigma(self):
        """
        Compute the covariance matrix of the MVN from the parameters

        """
        return flex.double([self.params[0] ** 2])

    def first_derivatives(self):
        """
        Compute the first derivatives of Sigma w.r.t the parameters

        """
        return flex.double([2 * self.params[0]])


class Angular2MosaicityParameterisation(BaseParameterisation):
    """
    A simple mosaicity parameterisation that uses 2 parameters to describe a
    multivariate normal angular mosaic spread. Sigma is enforced as positive
    definite by parameterising using the cholesky decomposition.
    W = | w1 0  0  |
        | 0 w1  0  |
        | 0  0 w2 |
    S = W*W^T
    """

    @staticmethod
    def is_angular():
        """
        Is angular

        """
        return True

    @staticmethod
    def num_parameters():
        """
        Get the number of parameters
        """
        return 2

    def sigma(self):
        """
        Compute the covariance matrix of the MVN from the parameters
        """
        M = matrix.sqr(
            (self.params[0], 0, 0, 0, self.params[0], 0, 0, 0, self.params[1])
        )
        return M * M.transpose()

    def first_derivatives(self):
        """
        Compute the first derivatives of Sigma w.r.t the parameters
        """
        b1, b2 = self.params

        d1 = (2 * b1, 0, 0, 0, 2 * b1, 0, 0, 0, 0)

        d2 = (0, 0, 0, 0, 0, 0, 0, 0, 2 * b2)

        return flex.mat3_double([d1, d2])


class Angular4MosaicityParameterisation(BaseParameterisation):
    """
    A simple mosaicity parameterisation that uses 4 parameters to describe a
    multivariate normal angular mosaic spread. Sigma is enforced as positive
    definite by parameterising using the cholesky decomposition.
    W = | w1  0  0  |
        | w2 w3  0  |
        | 0   0 w4 |
    S = W*W^T
    """

    @staticmethod
    def is_angular():
        """
        Is angular

        """
        return True

    @staticmethod
    def num_parameters():
        """
        Get the number of parameters
        """
        return 4

    def sigma(self):
        """
        Compute the covariance matrix of the MVN from the parameters
        """
        M = matrix.sqr(
            (
                self.params[0],
                0,
                0,
                self.params[1],
                self.params[2],
                0,
                0,
                0,
                self.params[3],
            )
        )
        return M * M.transpose()

    def first_derivatives(self):
        """
        Compute the first derivatives of Sigma w.r.t the parameters
        """
        b1, b2, b3, b4 = self.params

        d1 = (2 * b1, b2, 0, b2, 0, 0, 0, 0, 0)

        d2 = (0, b1, 0, b1, 2 * b2, 0, 0, 0, 0)

        d3 = (0, 0, 0, 0, 2 * b3, 0, 0, 0, 0)

        d4 = (0, 0, 0, 0, 0, 0, 0, 0, 2 * b4)

        return flex.mat3_double([d1, d2, d3, d4])


class ModelState(object):
    """
    A class to keep track of the model state

    """

    def __init__(
        self,
        experiment,
        mosaicity_parameterisation,
        wavelength_parameterisation=None,
        fix_mosaic_spread=False,
        fix_wavelength_spread=True,
        fix_unit_cell=False,
        fix_orientation=False,
    ):
        """
        Initialise the model state

        """

        # Save the crystal model
        self.experiment = experiment
        self.crystal = experiment.crystal

        # The U and P parameterisation
        self.U_parameterisation = CrystalOrientationParameterisation(self.crystal)
        self.B_parameterisation = CrystalUnitCellParameterisation(self.crystal)

        # The M and L parameterisations
        self.M_parameterisation = mosaicity_parameterisation
        self.L_parameterisation = wavelength_parameterisation

        # Set the flags to fix parameters
        self._is_mosaic_spread_fixed = fix_mosaic_spread
        self._is_wavelength_spread_fixed = fix_wavelength_spread
        self._is_unit_cell_fixed = fix_unit_cell
        self._is_orientation_fixed = fix_orientation

        # Check wavelength parameterisation
        if not self._is_wavelength_spread_fixed:
            assert self.L_parameterisation is not None

    def is_orientation_fixed(self):
        """
        Return whether the orientation is fixed

        """
        return self._is_orientation_fixed

    def is_unit_cell_fixed(self):
        """
        Return whether the unit cell is fixed

        """
        return self._is_unit_cell_fixed

    def is_mosaic_spread_fixed(self):
        """
        Return whether the mosaic spread is fixed

        """
        return self._is_mosaic_spread_fixed

    def is_mosaic_spread_angular(self):
        """
        Return whether the mosaic spread is angular

        """
        return self.M_parameterisation.is_angular()

    def is_wavelength_spread_fixed(self):
        """
        Return whether the wavelength spread is fixed

        """
        return self._is_wavelength_spread_fixed

    def get_unit_cell(self):
        """
        Get the crystal unit cell

        """
        return self.crystal.get_unit_cell()

    def get_U(self):
        """
        Get the crystal U matrix

        """
        return matrix.sqr(self.crystal.get_U())

    def get_B(self):
        """
        Get the crystal B matrix

        """
        return matrix.sqr(self.crystal.get_B())

    def get_A(self):
        """
        Get the crystal A matrix

        """
        return matrix.sqr(self.crystal.get_A())

    def get_M(self):
        """
        Get the Sigma M matrix

        """
        return self.M_parameterisation.sigma()

    def get_L(self):
        """
        Get the L sigma

        """
        if self.L_parameterisation is not None:
            return self.L_parameterisation.sigma()
        return flex.double()

    def get_U_params(self):
        """
        Get the U parameters

        """
        return flex.double(self.U_parameterisation.get_param_vals())

    def get_B_params(self):
        """
        Get the B parameters

        """
        return flex.double(self.B_parameterisation.get_param_vals())

    def get_M_params(self):
        """
        Get the M parameters

        """
        return flex.double(self.M_parameterisation.parameters)

    def get_L_params(self):
        """
        Get the L parameters

        """
        if self.L_parameterisation is not None:
            return flex.double(self.L_parameterisation.parameters)
        return []

    def set_U_params(self, params):
        """
        Set the U parameters

        """
        return self.U_parameterisation.set_param_vals(params)

    def set_B_params(self, params):
        """
        Set the B parameters

        """
        return self.B_parameterisation.set_param_vals(params)

    def set_M_params(self, params):
        """
        Set the M parameters

        """
        self.M_parameterisation.parameters = params

    def set_L_params(self, params):
        """
        Set the L parameters

        """
        if self.L_parameterisation is not None:
            self.L_parameterisation.parameters = params

    def num_U_params(self):
        """
        Get the number of U parameters

        """
        return len(self.get_U_params())

    def num_B_params(self):
        """
        Get the number of B parameters

        """
        return len(self.get_B_params())

    def num_M_params(self):
        """
        Get the number of M parameters

        """
        return len(self.get_M_params())

    def num_L_params(self):
        """
        Get the number of L parameters

        """
        return len(self.get_L_params())

    def get_dU_dp(self):
        """
        Get the first derivatives of U w.r.t its parameters

        """
        return flex.mat3_double(self.U_parameterisation.get_ds_dp())

    def get_dB_dp(self):
        """
        Get the first derivatives of B w.r.t its parameters

        """
        return flex.mat3_double(self.B_parameterisation.get_ds_dp())

    def get_dM_dp(self):
        """
        Get the first derivatives of M w.r.t its parameters

        """
        return self.M_parameterisation.first_derivatives()

    def get_dL_dp(self):
        """
        Get the first derivatives of L w.r.t its parameters

        """
        if self.L_parameterisation is not None:
            return self.L_parameterisation.first_derivatives()
        return flex.double()

    def get_active_parameters(self):
        """
        Get the active parameters in order: U, B, M, L, W

        """
        active_params = flex.double()
        if not self._is_orientation_fixed:
            active_params.extend(self.get_U_params())
        if not self._is_unit_cell_fixed:
            active_params.extend(self.get_B_params())
        if not self._is_mosaic_spread_fixed:
            active_params.extend(self.get_M_params())
        if not self._is_wavelength_spread_fixed:
            active_params.extend(self.get_L_params())
        assert len(active_params) > 0
        return active_params

    def set_active_parameters(self, params):
        """
        Set the active parameters in order: U, B, M, L, W

        """
        if not self._is_orientation_fixed:
            temp = params[: self.num_U_params()]
            params = params[self.num_U_params() :]
            self.set_U_params(temp)
        if not self._is_unit_cell_fixed:
            temp = params[: self.num_B_params()]
            params = params[self.num_B_params() :]
            self.set_B_params(temp)
        if not self._is_mosaic_spread_fixed:
            temp = params[: self.num_M_params()]
            params = params[self.num_M_params() :]
            self.set_M_params(temp)

        if not self._is_wavelength_spread_fixed:
            temp = params[: self.num_L_params()]
            params = params[self.num_L_params() :]
            self.set_L_params(temp)

    def get_labels(self):
        """
        Get the parameter labels

        """
        labels = []
        if not self._is_orientation_fixed:
            for i in range(len(self.get_U_params())):
                labels.append("Crystal_U_%d" % i)
        if not self._is_unit_cell_fixed:
            for i in range(len(self.get_B_params())):
                labels.append("Crystal_B_%d" % i)
        if not self._is_mosaic_spread_fixed:
            for i in range(len(self.get_M_params())):
                labels.append("Mosaicity_%d" % i)
        if not self._is_wavelength_spread_fixed:
            labels.append("Wavelength_Spread")
        assert len(labels) > 0
        return labels


class ReflectionModelState(object):
    """
    Class to compute basic derivatives of Sigma and r w.r.t parameters

    """

    def __init__(self, state, s0, h):
        """
        Initialise with the state and compute derivatives

        """
        # Get a load of matrices
        A = state.get_A()
        U = state.get_U()
        B = state.get_B()
        M = state.get_M()
        L = state.get_L()

        # Compute the reciprocal lattice vector
        h = matrix.col(h)
        r = A * h
        # t = r.length()

        # Define rotation for W sigma components
        q1 = r.cross(s0).normalize()
        q2 = r.cross(q1).normalize()
        q3 = r.normalize()
        Q = matrix.sqr(q1.elems + q2.elems + q3.elems)

        # rs0 = r.cross(s0)
        # rq1 = r.cross(q1)

        # Compute the covariance matrix
        if state.is_mosaic_spread_angular():
            self._sigma = Q.transpose() * M * Q
        else:
            self._sigma = M

        # Set the reciprocal lattice vector
        self._r = r

        # Get the derivatives of the model state
        dU_dp = state.get_dU_dp()
        dB_dp = state.get_dB_dp()
        dM_dp = state.get_dM_dp()
        dL_dp = state.get_dL_dp()

        # Get the wavelength spread
        assert len(L) == len(dL_dp)
        if len(L) == 0:
            self._sigma_lambda = 0
        else:
            assert len(L) == 1
            self._sigma_lambda = L[0]

        # The array of derivatives
        self._dr_dp = flex.vec3_double()
        self._ds_dp = flex.mat3_double()
        self._dl_dp = flex.double()

        # Compute derivatives w.r.t U parameters
        if not state.is_orientation_fixed():

            dr_dp_u = flex.vec3_double(state.num_U_params())
            ds_dp_u = flex.mat3_double(state.num_U_params())
            dl_dp_u = flex.double(state.num_U_params())
            for i in range(state.num_U_params()):
                dr_dp_u[i] = matrix.sqr(dU_dp[i]) * B * h
                # dr = matrix.col(dr_dp_u[i])
                # dt = r.dot(matrix.col(dr)) / t
                # drs0 = dr.cross(s0)
                # dq1 = drs0/rs0.length() - rs0*rs0.dot(drs0)/rs0.length()**3
                # drq1 = dr.cross(q1) + r.cross(dq1)
                # dq2 = drq1/rq1.length() - rq1*rq1.dot(drq1)/rq1.length()**3
                # dq3 = dr/r.length() - r*r.dot(dr)/r.length()**3
                # dQ = matrix.sqr(
                #   dq1.elems +
                #   dq2.elems +
                #   dq3.elems)
                # if state.is_mosaic_spread_angular():
                #   ds_dp_u[i] = Q.transpose()*M*Q \
                #              + dQ.transpose()*M*Q \
                #              + Q.transpose()*M*dQ
                # ds_dp_u[i] = dt*Q*(L+W)*Q.transpose() \
                #            + t*dQ*(L+W)*Q.transpose() \
                #            + t*Q*(L+W)*dQ.transpose()
            self._dr_dp.extend(dr_dp_u)
            self._ds_dp.extend(ds_dp_u)
            self._dl_dp.extend(dl_dp_u)

        # Compute derivatives w.r.t B parameters
        if not state.is_unit_cell_fixed():

            dr_dp_b = flex.vec3_double(state.num_B_params())
            ds_dp_b = flex.mat3_double(state.num_B_params())
            dl_dp_b = flex.double(state.num_B_params())
            for i in range(state.num_B_params()):
                dr_dp_b[i] = U * matrix.sqr(dB_dp[i]) * h
                # dr = matrix.col(dr_dp_b[i])
                # dt = r.dot(matrix.col(dr)) / t
                # drs0 = dr.cross(s0)
                # dq1 = drs0/rs0.length() - rs0*rs0.dot(drs0)/rs0.length()**3
                # drq1 = dr.cross(q1) + r.cross(dq1)
                # dq2 = drq1/rq1.length() - rq1*rq1.dot(drq1)/rq1.length()**3
                # dq3 = dr/r.length() - r*r.dot(dr)/r.length()**3
                # dQ = matrix.sqr(
                #   dq1.elems +
                #   dq2.elems +
                #   dq3.elems)
                # if state.is_mosaic_spread_angular():
                #   ds_dp_b[i] = Q.transpose()*M*Q \
                #              + dQ.transpose()*M*Q \
                # + Q.transpose()*M*dQ
                # ds_dp_b[i] = dt*Q*(L+W)*Q.transpose() \
                #            + t*dQ*(L+W)*Q.transpose() \
                #            + t*Q*(L+W)*dQ.transpose()
            self._dr_dp.extend(dr_dp_b)
            self._ds_dp.extend(ds_dp_b)
            self._dl_dp.extend(dl_dp_b)

        # Compute derivatives w.r.t M parameters
        if not state.is_mosaic_spread_fixed():
            dr_dp_m = flex.vec3_double(state.num_M_params())
            ds_dp_m = flex.mat3_double(state.num_M_params())
            dl_dp_m = flex.double(state.num_M_params())
            if state.is_mosaic_spread_angular():
                for i in range(state.num_M_params()):
                    dr_dp_m[i] = (0, 0, 0)
                    ds_dp_m[i] = Q.transpose() * matrix.sqr(dM_dp[i]) * Q
            else:
                for i in range(state.num_M_params()):
                    dr_dp_m[i] = (0, 0, 0)
                    ds_dp_m[i] = dM_dp[i]
            self._dr_dp.extend(dr_dp_m)
            self._ds_dp.extend(ds_dp_m)
            self._dl_dp.extend(dl_dp_m)

        # Compute derivatives w.r.t L parameters
        if not state.is_wavelength_spread_fixed():
            dr_dp_l = flex.vec3_double(state.num_L_params())
            ds_dp_l = flex.mat3_double(state.num_L_params())
            self._dr_dp.extend(dr_dp_l)
            self._ds_dp.extend(ds_dp_l)
            self._dl_dp.extend(dL_dp)

    def get_sigma(self):
        """
        Return the covariance matrix

        """
        return self._sigma

    def get_r(self):
        """
        Return the reciprocal lattice vector

        """
        return self._r

    def get_dS_dp(self):
        """
        Return the derivatives of the covariance matrix

        """
        return self._ds_dp

    def get_dr_dp(self):
        """
        Return the derivatives of the reciprocal lattice vector

        """
        return self._dr_dp

    def get_sigma_lambda(self):
        """
        Return the wavelength spread

        """
        return self._sigma_lambda

    def get_dL_dp(self):
        """
        Return the derivatives of the wavelength spread

        """
        return self._dl_dp
