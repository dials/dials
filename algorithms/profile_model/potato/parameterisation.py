from __future__ import division

from abc import ABC, abstractmethod
from typing import Dict, List, Optional

import numpy as np

from scitbx import linalg, matrix

from dials.algorithms.profile_model.potato import mosaicity_from_eigen_decomposition
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
            self.params = np.array([0.0] * self.num_parameters(), dtype=np.float)

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

    def sigma(self) -> matrix:
        """
        Compute the covariance matrix of the MVN from the parameters

        """
        M = matrix.sqr(
            (self.params[0], 0, 0, 0, self.params[0], 0, 0, 0, self.params[0])
        )
        return M * M.transpose()

    def first_derivatives(self) -> flex.mat3_double:
        """
        Compute the first derivatives of Sigma w.r.t the parameters

        """
        b1 = self.params[0]

        dSdb1 = (2 * b1, 0, 0, 0, 2 * b1, 0, 0, 0, 2 * b1)

        return flex.mat3_double([dSdb1])

    def mosaicity(self) -> Dict:
        """One value for mosaicity for Simple1"""
        m = mosaicity_from_eigen_decomposition(
            linalg.eigensystem.real_symmetric(
                self.sigma().as_flex_double_matrix()
            ).values()
        )
        return {"spherical": m[0]}


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

    def sigma(self) -> matrix:
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

    def mosaicity(self) -> Dict:
        """Three components for mosaicity for Simple6"""
        decomp = linalg.eigensystem.real_symmetric(self.sigma().as_flex_double_matrix())
        vals = list(mosaicity_from_eigen_decomposition(decomp.values()))
        min_m = min(vals)
        max_m = max(vals)
        vals.remove(min_m)
        vals.remove(max_m)
        return {"min": min_m, "mid": vals[0], "max": max_m}

    def first_derivatives(self) -> flex.mat3_double:
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

    def second_derivatives(self) -> flex.mat3_double:
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
        M = matrix.sqr(
            (self.params[0], 0, 0, 0, self.params[0], 0, 0, 0, self.params[1])
        )
        return M * M.transpose()

    def mosaicity(self) -> Dict:
        """Two unique components of mosaicity"""
        decomp = linalg.eigensystem.real_symmetric(self.sigma().as_flex_double_matrix())
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
    def is_angular() -> bool:
        return True

    @staticmethod
    def num_parameters() -> int:
        return 4

    def sigma(self) -> matrix:
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

    def mosaicity(self) -> Dict:
        """Three components of mosaicity"""
        decomp = linalg.eigensystem.real_symmetric(self.sigma().as_flex_double_matrix())
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

    def first_derivatives(self) -> flex.mat3_double:
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
    def U_matrix(self):
        return matrix.sqr(self.crystal.get_U())

    @property
    def B_matrix(self):
        return matrix.sqr(self.crystal.get_B())

    @property
    def A_matrix(self):
        return matrix.sqr(self.crystal.get_A())

    @property
    def mosaicity_covariance_matrix(self) -> matrix:
        return self._M_parameterisation.sigma()

    @property
    def wavelength_spread(self) -> flex.double:
        if self._L_parameterisation is not None:
            return self._L_parameterisation.sigma()
        return flex.double()

    @property
    def U_params(self) -> flex.double:
        """Get the parameters of the orientation parameterisation"""
        return flex.double(self._U_parameterisation.get_param_vals())

    @U_params.setter
    def U_params(self, params) -> None:
        self._U_parameterisation.set_param_vals(params)

    @property
    def B_params(self) -> flex.double:
        """Get the parameters of the orientation parameterisation"""
        return flex.double(self._B_parameterisation.get_param_vals())

    @B_params.setter
    def B_params(self, params) -> None:
        self._B_parameterisation.set_param_vals(params)

    @property
    def M_params(self) -> flex.double:
        "Parameters of the mosaicity parameterisation"
        return flex.double(self._M_parameterisation.parameters)

    @M_params.setter
    def M_params(self, params) -> None:
        self._M_parameterisation.parameters = params

    @property
    def L_params(self) -> flex.double:
        "Parameters of the Lambda (wavelength) parameterisation"
        if self._L_parameterisation is not None:
            return flex.double(self._L_parameterisation.parameters)
        return flex.double()

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
    def dM_dp(self) -> flex.mat3_double:
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
    def active_parameters(self) -> flex.double:
        """
        The active parameters in order: U, B, M, L, W
        """
        active_params = flex.double()
        if not self.is_orientation_fixed:
            active_params.extend(self.U_params)
        if not self.is_unit_cell_fixed:
            active_params.extend(self.B_params)
        if not self.is_mosaic_spread_fixed:
            active_params.extend(self.M_params)
        if not self.is_wavelength_spread_fixed:
            active_params.extend(self.L_params)
        assert len(active_params) > 0
        return active_params

    @active_parameters.setter
    def active_parameters(self, params) -> None:
        """
        Set the active parameters in order: U, B, M, L, W
        """
        if not self.is_orientation_fixed:
            n_U_params = self.U_params.size()
            temp = params[:n_U_params]
            params = params[n_U_params:]
            self.U_params = temp
        if not self.is_unit_cell_fixed:
            n_B_params = self.B_params.size()
            temp = params[:n_B_params]
            params = params[n_B_params:]
            self.B_params = temp
        if not self.is_mosaic_spread_fixed:
            n_M_params = self.M_params.size()
            temp = params[:n_M_params]
            params = params[n_M_params:]
            self.M_params = temp
        if not self.is_wavelength_spread_fixed:
            n_L_params = self.L_params.size()
            temp = params[:n_L_params]
            params = params[n_L_params:]
            self.L_params = temp

    @property
    def parameter_labels(self) -> List[str]:
        """
        Get the parameter labels

        """
        labels = []
        if not self.is_orientation_fixed:
            labels += [f"Crystal_U_{i}" for i in range(self.U_params.size())]
        if not self.is_unit_cell_fixed:
            labels += [f"Crystal_B_{i}" for i in range(self.B_params.size())]
        if not self.is_mosaic_spread_fixed:
            labels += [f"Mosaicity_{i}" for i in range(self.M_params.size())]
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
        h = matrix.col(h)
        r = state.A_matrix * h

        # Define rotation for W sigma components
        q1 = r.cross(s0).normalize()
        q2 = r.cross(q1).normalize()
        q3 = r.normalize()
        Q = matrix.sqr(q1.elems + q2.elems + q3.elems)
        # Compute the covariance matrix
        if state.is_mosaic_spread_angular:
            self._sigma = Q.transpose() * M * Q
        else:
            self._sigma = M

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
        self._dr_dp = flex.vec3_double()
        self._ds_dp = flex.mat3_double()
        self._dl_dp = flex.double()

        # Compute derivatives w.r.t U parameters
        if not state.is_orientation_fixed:
            n_U_params = state.U_params.size()
            dr_dp_u = flex.vec3_double(n_U_params)
            ds_dp_u = flex.mat3_double(n_U_params)
            dl_dp_u = flex.double(n_U_params)
            for i in range(n_U_params):
                dr_dp_u[i] = matrix.sqr(dU_dp[i]) * state.B_matrix * h
            self._dr_dp.extend(dr_dp_u)
            self._ds_dp.extend(ds_dp_u)
            self._dl_dp.extend(dl_dp_u)

        # Compute derivatives w.r.t B parameters
        if not state.is_unit_cell_fixed:
            n_B_params = state.B_params.size()
            dr_dp_b = flex.vec3_double(n_B_params)
            ds_dp_b = flex.mat3_double(n_B_params)
            dl_dp_b = flex.double(n_B_params)
            for i in range(n_B_params):
                dr_dp_b[i] = state.U_matrix * matrix.sqr(dB_dp[i]) * h
            self._dr_dp.extend(dr_dp_b)
            self._ds_dp.extend(ds_dp_b)
            self._dl_dp.extend(dl_dp_b)

        # Compute derivatives w.r.t M parameters
        if not state.is_mosaic_spread_fixed:
            n_M_params = state.M_params.size()
            dr_dp_m = flex.vec3_double(n_M_params)
            ds_dp_m = flex.mat3_double(n_M_params)
            dl_dp_m = flex.double(n_M_params)
            if state.is_mosaic_spread_angular:
                for i in range(n_M_params):
                    dr_dp_m[i] = (0, 0, 0)
                    ds_dp_m[i] = Q.transpose() * matrix.sqr(dM_dp[i]) * Q
            else:
                for i in range(n_M_params):
                    dr_dp_m[i] = (0, 0, 0)
                    ds_dp_m[i] = dM_dp[i]
            self._dr_dp.extend(dr_dp_m)
            self._ds_dp.extend(ds_dp_m)
            self._dl_dp.extend(dl_dp_m)

        # Compute derivatives w.r.t L parameters
        if not state.is_wavelength_spread_fixed:
            n_L_params = state.L_params.size()
            dr_dp_l = flex.vec3_double(n_L_params)
            ds_dp_l = flex.mat3_double(n_L_params)
            self._dr_dp.extend(dr_dp_l)
            self._ds_dp.extend(ds_dp_l)
            self._dl_dp.extend(dL_dp)

    @property
    def mosaicity_covariance_matrix(self) -> matrix:
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

    @property
    def wavelength_spread(self) -> float:
        return self._sigma_lambda

    def get_dL_dp(self):
        """
        Return the derivatives of the wavelength spread

        """
        return self._dl_dp
