from __future__ import annotations

from abc import ABC, abstractmethod
from typing import Dict, List, Optional

import numpy as np
from numpy.linalg import norm

from dxtbx import flumpy
from scitbx import linalg, matrix

from dials.algorithms.profile_model.ellipsoid import mosaicity_from_eigen_decomposition
from dials.algorithms.refinement.parameterisation.crystal_parameters import (
    CrystalOrientationParameterisation,
    CrystalUnitCellParameterisation,
)
from dials.array_family import flex


class BaseParameterisation(ABC):
    def __init__(self, params: Optional[np.array] = None) -> None:
        """
        Initialise with the parameters

        """
        if params is not None:
            assert len(params) == self.num_parameters()
            self.params = params
        else:
            self.params = np.array([0.0] * self.num_parameters(), dtype=np.float64)

    @abstractmethod
    def num_parameters(self):
        pass

    @property
    def parameters(self) -> np.array:
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
    def is_angular() -> bool:
        return False

    @staticmethod
    def num_parameters() -> int:
        return 1

    def sigma(self) -> np.array:
        """
        Compute the covariance matrix of the MVN from the parameters

        """
        M = np.array(
            [[self.params[0], 0, 0], [0, self.params[0], 0], [0, 0, self.params[0]]],
            dtype=np.float64,
        )
        return np.matmul(M, M.T)

    def first_derivatives(self) -> np.array:
        """
        Compute the first derivatives of Sigma w.r.t the parameters

        """
        b1 = self.params[0]

        # dSdb1 = (2 * b1, 0, 0, 0, 2 * b1, 0, 0, 0, 2 * b1)

        return np.array(
            [[[2.0 * b1, 0, 0], [0, 2.0 * b1, 0], [0, 0, 2.0 * b1]]], dtype=np.float64
        ).reshape(
            3, 3, 1
        )  # flex.mat3_double([dSdb1])

    def mosaicity(self) -> Dict:
        """One value for mosaicity for Simple1"""
        decomp = linalg.eigensystem.real_symmetric(
            matrix.sqr(flumpy.from_numpy(self.sigma())).as_flex_double_matrix()
        )
        v = list(mosaicity_from_eigen_decomposition(decomp.values()))
        return {"spherical": v[0]}


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
    def is_angular() -> bool:
        return False

    @staticmethod
    def num_parameters() -> int:
        return 6

    def sigma(self) -> np.array:
        """
        Compute the covariance matrix of the MVN from the parameters

        """
        M = np.array(
            [
                [self.params[0], 0, 0],
                [self.params[1], self.params[2], 0],
                [self.params[3], self.params[4], self.params[5]],
            ],
            dtype=np.float64,
        )
        return np.matmul(M, M.T)

    def mosaicity(self) -> Dict:
        """Three components for mosaicity for Simple6"""
        decomp = linalg.eigensystem.real_symmetric(
            matrix.sqr(flumpy.from_numpy(self.sigma())).as_flex_double_matrix()
        )
        vals = list(mosaicity_from_eigen_decomposition(decomp.values()))
        min_m = min(vals)
        max_m = max(vals)
        vals.remove(min_m)
        vals.remove(max_m)
        return {"min": min_m, "mid": vals[0], "max": max_m}

    def first_derivatives(self) -> np.array:
        """
        Compute the first derivatives of Sigma w.r.t the parameters

        """
        b1, b2, b3, b4, b5, b6 = self.params

        d1 = np.array(
            [[2 * b1, b2, b4], [b2, 0, 0], [b4, 0, 0]], dtype=np.float64
        ).reshape(3, 3, 1)

        d2 = np.array(
            [[0, b1, 0], [b1, 2 * b2, b4], [0, b4, 0]], dtype=np.float64
        ).reshape(3, 3, 1)

        d3 = np.array(
            [[0, 0, 0], [0, 2 * b3, b5], [0, b5, 0]], dtype=np.float64
        ).reshape(3, 3, 1)

        d4 = np.array(
            [[0, 0, b1], [0, 0, b2], [b1, b2, 2 * b4]], dtype=np.float64
        ).reshape(3, 3, 1)

        d5 = np.array(
            [[0, 0, 0], [0, 0, b3], [0, b3, 2 * b5]], dtype=np.float64
        ).reshape(3, 3, 1)

        d6 = np.array([[0, 0, 0], [0, 0, 0], [0, 0, 2 * b6]], dtype=np.float64).reshape(
            3, 3, 1
        )

        ds = np.concatenate([d1, d2, d3, d4, d5, d6], axis=2)

        return ds


class WavelengthSpreadParameterisation(BaseParameterisation):
    """
    A simple wavelength parameterisation that uses 1 parameter to describe a
    multivariate normal wavelength spread. Sigma is enforced as positive
    definite by parameterising using the cholesky decomposition.

    L = | 0 0 0  |
        | 0 0 0  |
        | 0 0 l1 |

    S = L*L^T

    """

    @staticmethod
    def num_parameters() -> int:
        return 1

    def sigma(self) -> flex.double:
        """
        The normal distribution sigma

        """
        return flex.double([self.params[0]])

    def first_derivatives(self) -> flex.double:
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
    def is_angular() -> bool:
        return True

    @staticmethod
    def num_parameters() -> int:
        return 2

    def sigma(self) -> matrix:
        """
        Compute the covariance matrix of the MVN from the parameters
        """
        M = np.array(
            [[self.params[0], 0, 0], [0, self.params[0], 0], [0, 0, self.params[1]]],
            dtype=np.float64,
        )
        return np.matmul(M, M.T)

    def mosaicity(self) -> Dict:
        """Two unique components of mosaicity"""
        decomp = linalg.eigensystem.real_symmetric(
            matrix.sqr(flumpy.from_numpy(self.sigma())).as_flex_double_matrix()
        )
        m = mosaicity_from_eigen_decomposition(decomp.values())
        v = decomp.vectors()
        mosaicities = {"radial": 0, "angular": 0}
        # two values must be same, could have accidental degeneracy where all 3 same:
        unique_ = list(set(m))
        if len(unique_) == 1:
            return {"angular": unique_[0], "radial": unique_[0]}
        else:
            assert len(unique_) == 2
            for i in range(3):
                vec = (v[i * 3], v[(i * 3) + 1], v[(i * 3) + 2])
                if vec == (0, 0, 1):
                    mosaicities["radial"] = m[i]
                else:
                    mosaicities["angular"] = m[i]
        return mosaicities

    def first_derivatives(self) -> flex.mat3_double:
        """
        Compute the first derivatives of Sigma w.r.t the parameters
        """
        b1, b2 = self.params

        d1 = np.array(
            [[2 * b1, 0, 0], [0, 2 * b1, 0], [0, 0, 0]], dtype=np.float64
        ).reshape(3, 3, 1)

        d2 = np.array([[0, 0, 0], [0, 0, 0], [0, 0, 2 * b2]], dtype=np.float64).reshape(
            3, 3, 1
        )

        return np.concatenate([d1, d2], axis=2)


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
    def is_angular() -> bool:
        return True

    @staticmethod
    def num_parameters() -> int:
        return 4

    def sigma(self) -> np.array:
        """
        Compute the covariance matrix of the MVN from the parameters
        """
        M = np.array(
            [
                [self.params[0], 0, 0],
                [self.params[1], self.params[2], 0],
                [0, 0, self.params[3]],
            ],
            dtype=np.float64,
        )
        return np.matmul(M, M.T)

    def mosaicity(self) -> Dict:
        """Three components of mosaicity"""
        decomp = linalg.eigensystem.real_symmetric(
            matrix.sqr(flumpy.from_numpy(self.sigma())).as_flex_double_matrix()
        )
        m = mosaicity_from_eigen_decomposition(decomp.values())
        v = decomp.vectors()
        mosaicities = {"radial": 0, "angular_0": 0, "angular_1": 0}
        n_angular = 0
        for i in range(3):
            vec = (v[i * 3], v[(i * 3) + 1], v[(i * 3) + 2])
            if vec == (0, 0, 1):
                mosaicities["radial"] = m[i]
            else:
                mosaicities["angular_" + str(n_angular)] = m[i]
                n_angular += 1
        return mosaicities

    def first_derivatives(self) -> np.array:
        """
        Compute the first derivatives of Sigma w.r.t the parameters
        """
        b1, b2, b3, b4 = self.params

        d1 = np.array(
            [[2 * b1, b2, 0], [b2, 0, 0], [0, 0, 0]], dtype=np.float64
        ).reshape(3, 3, 1)

        d2 = np.array(
            [[0, b1, 0], [b1, 2 * b2, 0], [0, 0, 0]], dtype=np.float64
        ).reshape(3, 3, 1)

        d3 = np.array([[0, 0, 0], [0, 2 * b3, 0], [0, 0, 0]], dtype=np.float64).reshape(
            3, 3, 1
        )

        d4 = np.array([[0, 0, 0], [0, 0, 0], [0, 0, 2 * b4]], dtype=np.float64).reshape(
            3, 3, 1
        )

        # return flex.mat3_double([d1, d2, d3, d4])

        return np.concatenate([d1, d2, d3, d4], axis=2)  # flex.mat3_double([dSdb1])


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
        self._U_parameterisation = CrystalOrientationParameterisation(self.crystal)
        self._B_parameterisation = CrystalUnitCellParameterisation(self.crystal)

        # The M and L parameterisations
        self._M_parameterisation = mosaicity_parameterisation
        self._L_parameterisation = wavelength_parameterisation

        # Set the flags to fix parameters
        self._is_mosaic_spread_fixed = fix_mosaic_spread
        self._is_wavelength_spread_fixed = fix_wavelength_spread
        self._is_unit_cell_fixed = fix_unit_cell
        self._is_orientation_fixed = fix_orientation

        # Check wavelength parameterisation
        if not self.is_wavelength_spread_fixed:
            assert self._L_parameterisation is not None

    @property
    def is_orientation_fixed(self) -> bool:
        return self._is_orientation_fixed

    @property
    def is_unit_cell_fixed(self) -> bool:
        return self._is_unit_cell_fixed

    @property
    def is_mosaic_spread_fixed(self) -> bool:
        return self._is_mosaic_spread_fixed

    @property
    def is_mosaic_spread_angular(self) -> bool:
        return self._M_parameterisation.is_angular()

    @property
    def is_wavelength_spread_fixed(self) -> bool:
        return self._is_wavelength_spread_fixed

    @property
    def unit_cell(self):
        return self.crystal.get_unit_cell()

    @property
    def U_matrix(self) -> np.array:
        return np.array([self.crystal.get_U()], dtype=np.float64).reshape(3, 3)

    @property
    def B_matrix(self) -> np.array:
        return np.array([self.crystal.get_B()], dtype=np.float64).reshape(3, 3)

    @property
    def A_matrix(self) -> np.array:
        return np.array([self.crystal.get_A()], dtype=np.float64).reshape(3, 3)

    @property
    def mosaicity_covariance_matrix(self) -> np.array:
        return self._M_parameterisation.sigma()

    @property
    def wavelength_spread(self) -> flex.double:
        if self._L_parameterisation is not None:
            return self._L_parameterisation.sigma()
        return flex.double()

    @property
    def U_params(self) -> np.array:
        """Get the parameters of the orientation parameterisation"""
        return np.array(self._U_parameterisation.get_param_vals(), dtype=np.float64)

    @U_params.setter
    def U_params(self, params) -> None:
        self._U_parameterisation.set_param_vals(tuple(float(i) for i in params))

    @property
    def B_params(self) -> np.array:
        """Get the parameters of the orientation parameterisation"""
        return np.array(self._B_parameterisation.get_param_vals(), dtype=np.float64)

    @B_params.setter
    def B_params(self, params) -> None:
        self._B_parameterisation.set_param_vals(tuple(float(i) for i in params))

    @property
    def M_params(self) -> np.array:
        "Parameters of the mosaicity parameterisation"
        return self._M_parameterisation.parameters

    @M_params.setter
    def M_params(self, params) -> None:
        self._M_parameterisation.parameters = params

    @property
    def L_params(self) -> np.array:
        "Parameters of the Lambda (wavelength) parameterisation"
        if self._L_parameterisation is not None:
            return np.array(self._L_parameterisation.parameters, dtype=np.float64)
        return np.array([])

    @L_params.setter
    def L_params(self, params: flex.double) -> None:
        if self._L_parameterisation is not None:
            self._L_parameterisation.parameters = params

    @property
    def dU_dp(self) -> flex.mat3_double:
        """
        Get the first derivatives of U w.r.t its parameters

        """
        return flex.mat3_double(self._U_parameterisation.get_ds_dp())

    @property
    def dB_dp(self) -> flex.mat3_double:
        """
        Get the first derivatives of B w.r.t its parameters

        """
        return flex.mat3_double(self._B_parameterisation.get_ds_dp())

    @property
    def dM_dp(self) -> np.array:
        """
        Get the first derivatives of M w.r.t its parameters

        """
        return self._M_parameterisation.first_derivatives()

    @property
    def dL_dp(self) -> flex.double:
        """
        Get the first derivatives of L w.r.t its parameters

        """
        if self._L_parameterisation is not None:
            return self._L_parameterisation.first_derivatives()
        return flex.double()

    @property
    def active_parameters(self) -> np.array:
        """
        The active parameters in order: U, B, M, L, W
        """
        active_params = []
        if not self.is_orientation_fixed:
            active_params.append(self.U_params)
        if not self.is_unit_cell_fixed:
            active_params.append(self.B_params)
        if not self.is_mosaic_spread_fixed:
            active_params.append(self.M_params)
        if not self.is_wavelength_spread_fixed:
            active_params.append(self.L_params)
        active_params = np.concatenate(active_params)
        assert len(active_params) > 0
        return active_params

    @active_parameters.setter
    def active_parameters(self, params) -> None:
        """
        Set the active parameters in order: U, B, M, L, W
        """
        if not self.is_orientation_fixed:
            n_U_params = len(self.U_params)
            temp = params[:n_U_params]
            params = params[n_U_params:]
            self.U_params = temp
        if not self.is_unit_cell_fixed:
            n_B_params = len(self.B_params)
            temp = params[:n_B_params]
            params = params[n_B_params:]
            self.B_params = temp
        if not self.is_mosaic_spread_fixed:
            n_M_params = len(self.M_params)
            temp = params[:n_M_params]
            params = params[n_M_params:]
            self.M_params = np.array(temp)
        if not self.is_wavelength_spread_fixed:
            n_L_params = self.L_params.size
            temp = params[:n_L_params]
            self.L_params = temp

    @property
    def parameter_labels(self) -> List[str]:
        """
        Get the parameter labels

        """
        labels = []
        if not self.is_orientation_fixed:
            labels += [f"Crystal_U_{i}" for i in range(self.U_params.size)]
        if not self.is_unit_cell_fixed:
            labels += [f"Crystal_B_{i}" for i in range(self.B_params.size)]
        if not self.is_mosaic_spread_fixed:
            labels += [f"Mosaicity_{i}" for i in range(self.M_params.size)]
        if not self.is_wavelength_spread_fixed:
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
        M = state.mosaicity_covariance_matrix
        L = state.wavelength_spread

        # Compute the reciprocal lattice vector
        h = np.array([h], dtype=np.float64).reshape(3, 1)
        r = np.matmul(state.A_matrix, h).reshape(3, 1)
        s0 = np.array([s0], dtype=np.float64).reshape(3, 1)

        # Define rotation for W sigma components
        q1 = np.cross(r.flatten(), s0.flatten())
        q1 /= norm(q1)
        q2 = np.cross(r.flatten(), q1)
        q2 /= norm(q2)
        q3 = r / norm(r)
        q3 = q3.flatten()
        Q = np.array([q1, q2, q3], dtype=np.float64).reshape(3, 3)

        # Compute the covariance matrix
        if state.is_mosaic_spread_angular:
            self._sigma = np.matmul(np.matmul(Q.T, M), Q)
        else:
            self._sigma = M  #

        # Set the reciprocal lattice vector
        self._r = r

        # Get the derivatives of the model state
        dU_dp = state.dU_dp
        dB_dp = state.dB_dp
        dM_dp = state.dM_dp
        dL_dp = state.dL_dp

        # Get the wavelength spread
        assert len(L) == len(dL_dp)
        if len(L) == 0:
            self._sigma_lambda = 0
        else:
            assert len(L) == 1
            self._sigma_lambda = L[0]

        # The array of derivatives
        self._dr_dp = np.array([], dtype=np.float64).reshape(3, 0)
        self._ds_dp = np.array([], dtype=np.float64).reshape(3, 3, 0)
        self._dl_dp = flex.double()

        # Compute derivatives w.r.t U parameters
        if not state.is_orientation_fixed:
            n_U_params = len(state.U_params)
            dr_dp_u = np.zeros(shape=(3, n_U_params), dtype=np.float64)
            ds_dp_u = np.zeros(shape=(3, 3, n_U_params), dtype=np.float64)
            dl_dp_u = flex.double(n_U_params)
            for i in range(n_U_params):
                v = np.matmul(
                    np.matmul(
                        np.array(dU_dp[i], dtype=np.float64).reshape(3, 3),
                        state.B_matrix,
                    ),
                    h,
                )
                dr_dp_u[:, i] = v[:, 0]
            self._dr_dp = np.concatenate([self._dr_dp, dr_dp_u], axis=1)
            self._ds_dp = np.concatenate([self._ds_dp, ds_dp_u], axis=2)
            self._dl_dp.extend(dl_dp_u)

        # Compute derivatives w.r.t B parameters
        if not state.is_unit_cell_fixed:
            n_B_params = len(state.B_params)
            dr_dp_b = np.zeros(shape=(3, n_B_params), dtype=np.float64)
            ds_dp_b = np.zeros(shape=(3, 3, n_B_params), dtype=np.float64)
            dl_dp_b = flex.double(n_B_params)
            for i in range(n_B_params):
                dr_dp_b[:, i] = np.matmul(
                    np.matmul(
                        state.U_matrix,
                        np.array(dB_dp[i], dtype=np.float64).reshape(3, 3),
                    ),
                    h,
                )[:, 0]
            self._dr_dp = np.concatenate([self._dr_dp, dr_dp_b], axis=1)
            self._ds_dp = np.concatenate([self._ds_dp, ds_dp_b], axis=2)
            self._dl_dp.extend(dl_dp_b)

        # Compute derivatives w.r.t M parameters
        if not state.is_mosaic_spread_fixed:
            n_M_params = len(state.M_params)
            dr_dp_m = np.zeros(shape=(3, n_M_params), dtype=np.float64)
            ds_dp_m = np.zeros(shape=(3, 3, n_M_params), dtype=np.float64)
            dl_dp_m = flex.double(n_M_params)
            if state.is_mosaic_spread_angular:
                for i in range(n_M_params):
                    ds_dp_m[:, :, i] = np.matmul(
                        np.matmul(
                            Q.T,
                            dM_dp[:, :, i],
                        ),
                        Q,
                    )
            else:
                for i in range(n_M_params):
                    ds_dp_m[:, :, i] = dM_dp[:, :, i]
            self._dr_dp = np.concatenate([self._dr_dp, dr_dp_m], axis=1)
            self._ds_dp = np.concatenate([self._ds_dp, ds_dp_m], axis=2)
            self._dl_dp.extend(dl_dp_m)
        # Compute derivatives   w.r.t L parameters
        if not state.is_wavelength_spread_fixed:
            n_l_params = len(state.L_params)
            dr_dp_l = np.zeros(shape=(3, n_l_params), dtype=np.float64)
            ds_dp_l = np.zeros(shape=(3, 3, n_l_params), dtype=np.float64)
            self._dr_dp = np.concatenate([self._dr_dp, dr_dp_l], axis=1)
            self._ds_dp = np.concatenate([self._ds_dp, ds_dp_l], axis=2)
            self._dl_dp.extend(dL_dp)

    @property
    def mosaicity_covariance_matrix(self) -> np.array:
        return self._sigma

    def get_r(self) -> np.array:
        """
        Return the reciprocal lattice vector

        """
        return self._r

    def get_dS_dp(self) -> np.array:
        """
        Return the derivatives of the covariance matrix (an array of size 3x3xn)

        """
        return self._ds_dp

    def get_dr_dp(self) -> np.array:
        """
        Return the derivatives of the reciprocal lattice vector ( an array of size 3xn)

        """
        return self._dr_dp

    @property
    def wavelength_spread(self) -> float:
        return self._sigma_lambda

    def get_dL_dp(self) -> flex.double:
        """
        Return the derivatives of the wavelength spread

        """
        return self._dl_dp
