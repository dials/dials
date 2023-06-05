from __future__ import annotations

from math import exp, sqrt

import numpy as np
from numpy.linalg import norm

from dxtbx import flumpy
from libtbx.phil import parse
from scitbx import linalg, matrix
from scitbx.linalg import eigensystem, l_l_transpose_cholesky_decomposition_in_place

from dials.algorithms.profile_model.ellipsoid import (
    BBoxCalculatorAngular,
    BBoxCalculatorSimple,
    MaskCalculatorAngular,
    MaskCalculatorSimple,
    PredictorAngular,
    PredictorSimple,
)
from dials.algorithms.profile_model.ellipsoid.parameterisation import (
    Angular2MosaicityParameterisation,
    Angular4MosaicityParameterisation,
    Simple1Angular1MosaicityParameterisation,
    Simple1MosaicityParameterisation,
    Simple6Angular1MosaicityParameterisation,
    Simple6Angular3MosaicityParameterisation,
    Simple6MosaicityParameterisation,
)
from dials.array_family import flex
from dials.constants import FULL_PARTIALITY
from dials.model.experiment.profile import ProfileModelExt

phil_scope = parse(
    """
rlp_mosaicity {

    model = simple1 simple6 *simple1angular1 simple6angular1 simple6angular3
    .type = choice

}

wavelength_spread {

    model = *delta
    .type = choice

}

unit_cell {

    fixed = False
    .type = bool

}

orientation {

    fixed = False
    .type = bool

}

indexing {

    fail_on_bad_index = False
      .type = bool

  }

refinement {

    max_separation = 2
        .type = float

    outlier_probability = 0.975
        .type = float

    n_macro_cycles = 3
        .type = int

    n_cycles = 3
        .type = int

    min_n_reflections=10
        .type = int

    max_iter=1000
        .type = int
        .help = "Max number of iterations per refinement cycle"

    LL_tolerance=1e-6
        .type = float
        .help = "Convergence tolerance for log likelihood during refinement"

}

prediction {
    d_min = None
        .type = float

    probability = %f
        .type = float
}"""
    % FULL_PARTIALITY
)


class EllipsoidProfileModel(ProfileModelExt):

    """
    An overall model class that conforms to the requirements of a
    dxtbx.profile_model entry point.
    """

    name = "ellipsoid"

    def __init__(self, parameterisation):
        self.parameterisation = parameterisation

    @classmethod
    def create(
        cls, params, reflections, crystal, beam, detector, goniometer=None, scan=None
    ):
        # a method to allow model creation for the standard integration program
        # need to work out the sigma d i.e. do initial integration
        from dials.algorithms.profile_model.gaussian_rs.calculator import (
            ComputeEsdBeamDivergence,
        )

        if not reflections:
            raise ValueError(
                "Reflections needed to determine sigma_d for the ellipsoid integrator"
            )

        sel = reflections.get_flags(reflections.flags.strong)
        strong_refls = reflections.select(sel)
        # Compute and initial spot size estimate and beam vector
        sigma_d = ComputeEsdBeamDivergence(detector, strong_refls).sigma()
        return cls.from_sigma_d(params.ellipsoid.rlp_mosaicity.model, sigma_d)

    def compute_bbox(
        self, reflections, crystal, beam, detector, goniometer=None, scan=None, **kwargs
    ):
        raise ValueError(
            "Ellipsoid profile modelling not implemented outside of dials.ssx_integrate"
        )

    @classmethod
    def from_sigma_d(cls, model, sigma_d):
        if model == "simple1":
            return cls(Simple1ProfileModel.from_sigma_d(sigma_d))
        elif model == "simple6":
            return cls(Simple6ProfileModel.from_sigma_d(sigma_d))
        elif model == "simple1angular1":
            return cls(Simple1Angular1ProfileModel.from_sigma_d(sigma_d))
        elif model == "simple6angular1":
            return cls(Simple6Angular1ProfileModel.from_sigma_d(sigma_d))
        elif model == "simple6angular3":
            return cls(Simple6Angular3ProfileModel.from_sigma_d(sigma_d))
        elif model == "angular2":
            return cls(Angular2ProfileModel.from_sigma_d(sigma_d))
        elif model == "angular4":
            return cls(Angular4ProfileModel.from_sigma_d(sigma_d))

        raise RuntimeError(f"Unknown profile model: {model}")

    @classmethod
    def from_dict(cls, d):
        if d["parameterisation"] == "simple1":
            return cls(Simple1ProfileModel.from_dict(d))
        if d["parameterisation"] == "simple6":
            return cls(Simple6ProfileModel.from_dict(d))
        if d["parameterisation"] == "simple1angular1":
            return cls(Simple1Angular1ProfileModel.from_dict(d))
        if d["parameterisation"] == "simple6angular1":
            return cls(Simple6Angular1ProfileModel.from_dict(d))
        if d["parameterisation"] == "simple6angular3":
            return cls(Simple6Angular3ProfileModel.from_dict(d))
        if d["parameterisation"] == "angular2":
            return cls(Angular2ProfileModel.from_dict(d))
        if d["parameterisation"] == "angular4":
            return cls(Angular4ProfileModel.from_dict(d))
        raise RuntimeError(
            f"Unknown profile model parameterisation: {d['parameterisation']}"
        )

    def mosaicity(self):
        return self.parameterisation.mosaicity()

    def sigma(self):
        return self.parameterisation.sigma()

    def to_dict(self):
        d = self.parameterisation.to_dict()
        d["parameterisation"] = d.pop("__id__")
        d["__id__"] = "ellipsoid"
        return d


class ProfileModelBase(object):
    """
    Class to store profile model

    """

    def __init__(self, params):
        """
        Initialise the class

        """
        self.params = params
        self._n_obs = None

    def sigma(self):
        """
        Get the sigma

        """
        return self.parameterisation().sigma()

    @property
    def n_obs(self):
        return self._n_obs

    @n_obs.setter
    def n_obs(self, n_obs_):
        self._n_obs = n_obs_

    def update_model_state_parameters(self, state):
        """
        Update the model state with the parameters

        """
        state.set_M_params(self.params)

    def update_model(self, state):
        """
        Update the model

        """

        # Compute the eigen decomposition of the covariance matrix and check
        # largest eigen value
        sqr_mat = matrix.sqr(flumpy.from_numpy(self.sigma()))
        eigen_decomposition = eigensystem.real_symmetric(
            sqr_mat.as_flex_double_matrix()
        )
        L = eigen_decomposition.values()
        if L[0] > 1e-5:
            raise RuntimeError("Mosaicity matrix is unphysically large")

        self.params = state.M_params

    def to_dict(self):
        """Convert the model to a dictionary."""
        params = list(self.parameterisation().parameters)
        sigma = self.sigma()
        return {
            "__id__": self.__class__.name,
            "parameters": params,
            "sigma": sigma.tolist(),
            "n_obs": self._n_obs,
        }

    @classmethod
    def from_dict(cls, d):
        """Convert the model to a dictionary."""
        model = cls(d["parameters"])
        return model

    def mosaicity(self):
        return self.parameterisation().mosaicity()


class SimpleProfileModelBase(ProfileModelBase):
    """
    Base class for simple profile models

    """

    def predict_reflections(
        self, experiments, miller_indices, probability=FULL_PARTIALITY
    ):
        """
        Predict the reflections

        """
        predictor = PredictorSimple(
            experiments[0], matrix.sqr(flumpy.from_numpy(self.sigma())), probability
        )
        return predictor.predict(miller_indices)

    def compute_bbox(self, experiments, reflections, probability=FULL_PARTIALITY):
        """
        Compute the bounding box

        """
        calculator = BBoxCalculatorSimple(
            experiments[0], matrix.sqr(flumpy.from_numpy(self.sigma())), probability, 4
        )
        calculator.compute(reflections)

    def compute_mask(self, experiments, reflections, probability=FULL_PARTIALITY):
        """
        Compute the mask

        """
        calculator = MaskCalculatorSimple(
            experiments[0], matrix.sqr(flumpy.from_numpy(self.sigma())), probability
        )
        calculator.compute(reflections)

    def sigma_for_reflection(self, s0, r):
        """
        Get sigma for a reflections

        """
        return np.array(self.sigma()).reshape(3, 3)

    def compute_partiality(self, experiments, reflections):
        """
        Compute the partiality

        """
        s0 = np.array([experiments[0].beam.get_s0()], dtype=np.float64).reshape(3, 1)
        s0_length = norm(s0)
        n_obs = experiments[0].crystal.mosaicity.parameterisation.n_obs
        assert n_obs is not None
        ## partiality is defined relative to max possible observation, which is the plane
        # through the centre of the RLP perpendicular to the min variance
        # need the min variance, so do decomposition
        eigen_decomposition = linalg.eigensystem.real_symmetric(
            matrix.sqr(
                experiments[0].crystal.mosaicity.sigma().flatten()
            ).as_flex_double_matrix()
        )
        eigen_values = eigen_decomposition.values()
        S00 = min(eigen_values)
        partiality = flex.double(reflections.size())
        partiality_variance = flex.double(reflections.size())
        for k, s2_vec in enumerate(reflections["s2"]):
            s2 = np.array(list(s2_vec), dtype=np.float64).reshape(3, 1)
            sigma = experiments[0].crystal.mosaicity.sigma()
            R = compute_change_of_basis_operation(s0, s2)

            S = np.matmul(R, np.array(sigma).reshape(3, 3))
            S = np.matmul(S, R.T)
            mu = np.matmul(R, s2)

            mu_norm = mu / norm(mu)
            assert abs(1.0 - mu_norm.flatten()[2]) < 1e-7
            S22 = S[2, 2]
            mu2 = mu.flatten()[2]
            eps = s0_length - mu2
            var_eps = S22 / n_obs  # Approximation
            partiality[k] = exp(-0.5 * eps * (1 / S22) * eps) * sqrt(S00 / S22)
            partiality_variance[k] = (
                var_eps * (eps**2 / (S00 * S22)) * exp(eps**2 / S22)
            )

        reflections["partiality"] = partiality
        reflections["partiality.inv.variance"] = partiality_variance

    @classmethod
    def from_params(Class, params):
        """
        Create the class from some parameters

        """
        return Class(params)


class Simple1ProfileModel(SimpleProfileModelBase):
    """
    Simple 1 profile model class

    """

    name = "simple1"

    def parameterisation(self):
        """
        Get the parameterisation

        """
        return Simple1MosaicityParameterisation(self.params)

    @classmethod
    def from_sigma_d(Class, sigma_d):
        """
        Create the profile model from sigma_d estimate

        """
        return Class.from_params(np.array([sigma_d], dtype=np.float64))

    @classmethod
    def from_sigma(Class, sigma):
        """
        Construct the profile model from the sigma

        """

        # Construct triangular matrix
        LL = flex.double()
        for j in range(3):
            for i in range(j + 1):
                LL.append(sigma[j * 3 + i])

        # Do the cholesky decomposition
        _ = l_l_transpose_cholesky_decomposition_in_place(LL)
        TINY = 1e-6  ###FIXME
        assert abs(LL[1] - 0) < TINY
        assert abs(LL[2] - LL[0]) < TINY
        assert abs(LL[3] - 0) < TINY
        assert abs(LL[4] - 0) < TINY
        assert abs(LL[5] - LL[0]) < TINY

        # Setup the parameters
        return Class.from_params(flex.double((LL[0],)))


class Simple6ProfileModel(SimpleProfileModelBase):
    """
    Class to store profile model

    """

    name = "simple6"

    def parameterisation(self):
        """
        Get the parameterisation

        """
        return Simple6MosaicityParameterisation(self.params)

    @classmethod
    def from_sigma_d(Class, sigma_d):
        """
        Create the profile model from sigma_d estimate

        """
        return Class.from_params(
            np.array([sigma_d, 0, sigma_d, 0, 0, sigma_d], dtype=np.float64)
        )

    @classmethod
    def from_sigma(Class, sigma):
        """
        Construct the profile model from the sigma

        """

        # Construct triangular matrix
        LL = flex.double()
        for j in range(3):
            for i in range(j + 1):
                LL.append(sigma[j * 3 + i])

        # Do the cholesky decomposition
        _ = l_l_transpose_cholesky_decomposition_in_place(LL)

        # Setup the parameters
        return Class.from_params(
            flex.double((LL[0], LL[1], LL[2], LL[3], LL[4], LL[5]))
        )


class AngularProfileModelBase(ProfileModelBase):
    """
    Class to store profile model

    """

    def sigma_for_reflection(self, s0, r):
        """
        Sigma for a reflection

        """
        Q = compute_change_of_basis_operation(s0, r)
        scale = norm(r) ** 2
        S = np.array([[scale, 0.0, 0.0, 0.0, scale, 0.0, 0.0, 0.0, 0.0]]).reshape(3, 3)
        sigma = (
            np.matmul(
                np.matmul(
                    Q.T, S * np.array(self.parameterisation().sigma_A()).reshape(3, 3)
                ),
                Q,
            )
            + self.sigma()
        )
        return sigma

    def predict_reflections(
        self, experiments, miller_indices, probability=FULL_PARTIALITY
    ):
        """
        Predict the reflections

        """
        predictor = PredictorAngular(
            experiments[0],
            matrix.sqr(flumpy.from_numpy(self.sigma())),
            probability,
            matrix.sqr(flumpy.from_numpy(self.parameterisation().sigma_A())),
        )
        return predictor.predict(miller_indices)

    def compute_bbox(self, experiments, reflections, probability=FULL_PARTIALITY):
        """
        Compute the bounding box

        """
        calculator = BBoxCalculatorAngular(
            experiments[0],
            matrix.sqr(flumpy.from_numpy(self.sigma())),
            probability,
            4,
            matrix.sqr(flumpy.from_numpy(self.parameterisation().sigma_A())),
        )
        calculator.compute(reflections)

    def compute_mask(self, experiments, reflections, probability=FULL_PARTIALITY):
        """
        Compute the mask

        """
        calculator = MaskCalculatorAngular(
            experiments[0],
            matrix.sqr(flumpy.from_numpy(self.sigma())),
            probability,
            matrix.sqr(flumpy.from_numpy(self.parameterisation().sigma_A())),
        )
        calculator.compute(reflections)

    def compute_partiality(self, experiments, reflections):
        """
        Compute the partiality

        """
        """s0 = np.array([experiments[0].beam.get_s0()], dtype=np.float64).reshape(3, 1)
        s0_length = norm(s0)
        n_obs = experiments[0].crystal.mosaicity.parameterisation.n_obs
        assert n_obs is not None"""
        """sigma = experiments[0].crystal.mosaicity.sigma()
        sigma = np.array(sigma).reshape(3, 3)
        x, y, z = reflections["s2"].parts()
        s2 = np.array([x, y, z])
        ## partiality is defined relative to max possible observation, which is the plane
        # through the centre of the RLP perpendicular to the min variance
        # need the min variance, so do decomposition
        eigen_decomposition = linalg.eigensystem.real_symmetric(
            matrix.sqr(sigma.flatten()).as_flex_double_matrix()
        )
        eigen_values = eigen_decomposition.values()
        S00 = min(eigen_values)
        r = s2 - s0
        Rs = compute_change_of_basis_operations(s0, s2)
        Qs = compute_change_of_basis_operations(s0, r)
        sigma_qs = np.einsum("mda,db,mbc->mac", Qs, sigma, Qs)  # Q.T * sigma * Q
        Ss = np.einsum("mil,mlj,mkj->mik", Rs, sigma_qs, Rs)  # R * sigma_q * R.T
        mus = np.einsum("mij,jm->mi", Rs, s2)
        eps = mus[:, 2] - s0_length
        eps2 = np.square(eps)
        S22 = Ss[:, 2, 2]

        partiality = np.exp(-0.5 * eps2 / S22) * np.sqrt(S00 / S22)
        partiality_variance = eps2 * np.exp(eps2 / S22) / (n_obs * S00)

        reflections["partiality"] = flumpy.from_numpy(partiality)
        reflections["partiality.inv.variance"] = flumpy.from_numpy(partiality_variance)"""

        s0 = np.array([experiments[0].beam.get_s0()], dtype=np.float64).reshape(3, 1)
        s0_length = norm(s0)
        n_obs = experiments[0].crystal.mosaicity.parameterisation.n_obs
        assert n_obs is not None
        ## partiality is defined relative to max possible observation, which is the plane
        # through the centre of the RLP perpendicular to the min variance at r==0
        # need the min variance, so do decomposition
        eigen_decomposition = linalg.eigensystem.real_symmetric(
            matrix.sqr(self.sigma().flatten()).as_flex_double_matrix()
        )
        eigen_values = eigen_decomposition.values()
        S00 = min(eigen_values)
        partiality = flex.double(reflections.size())
        partiality_variance = flex.double(reflections.size())
        for k, s2_vec in enumerate(reflections["s2"]):
            s2 = np.array(list(s2_vec), dtype=np.float64).reshape(3, 1)
            r = s2 - s0
            sigma = self.sigma_for_reflection(
                s0, r
            )  #  experiments[0].crystal.mosaicity.sigma()
            R = compute_change_of_basis_operation(s0, s2)

            S = np.matmul(R, np.array(sigma).reshape(3, 3))
            S = np.matmul(S, R.T)
            mu = np.matmul(R, s2)

            mu_norm = mu / norm(mu)
            assert abs(1.0 - mu_norm.flatten()[2]) < 1e-7
            S22 = S[2, 2]
            mu2 = mu.flatten()[2]
            eps = s0_length - mu2
            var_eps = S22 / n_obs  # Approximation
            partiality[k] = exp(-0.5 * eps * (1 / S22) * eps) * sqrt(S00 / S22)
            partiality_variance[k] = (
                var_eps * (eps**2 / (S00 * S22)) * exp(eps**2 / S22)
            )

        reflections["partiality"] = partiality
        reflections["partiality.inv.variance"] = partiality_variance

    @classmethod
    def from_params(Class, params):
        """
        Create the class from some parameters

        """
        return Class(params)


class Simple1Angular1ProfileModel(AngularProfileModelBase):

    name = "simple1angular1"

    def parameterisation(self):
        """
        Get the parameterisation

        """
        return Simple1Angular1MosaicityParameterisation(self.params)

    @classmethod
    def from_sigma_d(Class, sigma_d):
        """
        Create the profile model from sigma_d estimate

        """
        return Class.from_params(np.array([sigma_d, sigma_d], dtype=np.float64))

    @classmethod
    def from_sigma(Class, sigma):
        """
        Construct the profile model from the sigma

        """

        # Construct triangular matrix
        LL = flex.double()
        for j in range(3):
            for i in range(j + 1):
                LL.append(sigma[j * 3 + i])

        # Do the cholesky decomposition
        _ = l_l_transpose_cholesky_decomposition_in_place(LL)

        # Check the sigma is as we expect
        TINY = 1e-10
        assert abs(LL[1] - 0) < TINY
        assert abs(LL[2] - LL[0]) < TINY
        assert abs(LL[3] - 0) < TINY
        assert abs(LL[4] - 0) < TINY

        # Setup the parameters
        return Class.from_params(flex.double((LL[0], LL[5])))


class Simple6Angular1ProfileModel(AngularProfileModelBase):

    name = "simple6angular1"

    def parameterisation(self):
        """
        Get the parameterisation

        """
        return Simple6Angular1MosaicityParameterisation(self.params)

    @classmethod
    def from_sigma_d(Class, sigma_d):
        """
        Create the profile model from sigma_d estimate

        """
        return Class.from_params(
            np.array(
                [sigma_d, 0, sigma_d, 0, 0, sigma_d, sigma_d / 10.0], dtype=np.float64
            )
        )

    @classmethod
    def from_sigma(Class, sigma):
        """
        Construct the profile model from the sigma

        """

        # Construct triangular matrix
        LL = flex.double()
        for j in range(3):
            for i in range(j + 1):
                LL.append(sigma[j * 3 + i])

        # Do the cholesky decomposition
        _ = l_l_transpose_cholesky_decomposition_in_place(LL)

        # Setup the parameters
        return Class.from_params(
            flex.double((LL[0], LL[1], LL[2], LL[3], LL[4], LL[5], LL[5]))
        )


class Simple6Angular3ProfileModel(AngularProfileModelBase):

    name = "simple6angular3"

    def parameterisation(self):
        """
        Get the parameterisation

        """
        return Simple6Angular3MosaicityParameterisation(self.params)

    @classmethod
    def from_sigma_d(Class, sigma_d):
        """
        Create the profile model from sigma_d estimate

        """
        return Class.from_params(
            np.array(
                [
                    sigma_d,
                    0,
                    sigma_d,
                    0,
                    0,
                    sigma_d,
                    sigma_d / 10.0,
                    0.0,
                    sigma_d / 10.0,
                ],
                dtype=np.float64,
            )
        )

    @classmethod
    def from_sigma(Class, sigma):
        """
        Construct the profile model from the sigma

        """

        # Construct triangular matrix
        LL = flex.double()
        for j in range(3):
            for i in range(j + 1):
                LL.append(sigma[j * 3 + i])

        # Do the cholesky decomposition
        _ = l_l_transpose_cholesky_decomposition_in_place(LL)

        # Setup the parameters
        return Class.from_params(
            flex.double((LL[0], LL[1], LL[2], LL[3], LL[4], LL[5], LL[5], 0, LL[5]))
        )


class Angular2ProfileModel(AngularProfileModelBase):
    """
    Class to store profile model

    """

    name = "angular2"

    def parameterisation(self):
        """
        Get the parameterisation

        """
        return Angular2MosaicityParameterisation(self.params)

    @classmethod
    def from_sigma_d(Class, sigma_d):
        """
        Create the profile model from sigma_d estimate

        """
        return Class.from_params(np.array([sigma_d, sigma_d], dtype=np.float64))

    @classmethod
    def from_sigma(Class, sigma):
        """
        Construct the profile model from the sigma

        """

        # Construct triangular matrix
        LL = flex.double()
        for j in range(3):
            for i in range(j + 1):
                LL.append(sigma[j * 3 + i])

        # Do the cholesky decomposition
        _ = l_l_transpose_cholesky_decomposition_in_place(LL)

        # Check the sigma is as we expect
        TINY = 1e-10
        assert abs(LL[1] - 0) < TINY
        assert abs(LL[2] - LL[0]) < TINY
        assert abs(LL[3] - 0) < TINY
        assert abs(LL[4] - 0) < TINY

        # Setup the parameters
        return Class.from_params(flex.double((LL[0], LL[5])))


class Angular4ProfileModel(AngularProfileModelBase):
    """
    Class to store profile model

    """

    name = "angular4"

    def parameterisation(self):
        """
        Get the parameterisation

        """
        return Angular4MosaicityParameterisation(self.params)

    @classmethod
    def from_sigma_d(Class, sigma_d):
        """
        Create the profile model from sigma_d estimate

        """
        return Class.from_params(
            np.array([sigma_d, 0, sigma_d, sigma_d], dtype=np.float64)
        )

    @classmethod
    def from_sigma(Class, sigma):
        """
        Construct the profile model from the sigma

        """

        # Construct triangular matrix
        LL = flex.double()
        for j in range(3):
            for i in range(j + 1):
                LL.append(sigma[j * 3 + i])

        # Do the cholesky decomposition
        _ = l_l_transpose_cholesky_decomposition_in_place(LL)

        # Check the sigma is as we expect
        TINY = 1e-10
        assert abs(LL[3] - 0) < TINY
        assert abs(LL[4] - 0) < TINY

        # Setup the parameters
        return Class.from_params(flex.double((LL[0], LL[1], LL[2], LL[5])))


class ProfileModelFactory(object):
    """
    Class to create profile models

    """

    @classmethod
    def from_sigma_d(Class, model, sigma_d):
        """
        Construct a profile model from an initial sigma estimate

        """
        return EllipsoidProfileModel.from_sigma_d(model, sigma_d)


def compute_change_of_basis_operations(s0, s2_array):
    s0 = s0.reshape(1, 3)
    s2 = s2_array.T
    assert s2.shape[1] == 3
    e1 = np.cross(s2, s0)
    e2 = np.cross(s2, e1)
    norm_e1 = norm(e1, axis=1).reshape(-1, 1)
    norm_e2 = norm(e2, axis=1).reshape(-1, 1)
    e1 /= norm_e1
    e2 /= norm_e2
    norm_s2 = norm(s2, axis=1).reshape(-1, 1)
    e3 = s2 / norm_s2
    n = e1.shape[0]
    R = np.hstack([e1, e2, e3]).reshape(n, 3, 3)
    return R


def compute_change_of_basis_operation(s0, s2):
    s2 = s2.flatten()
    s0 = s0.flatten()
    e1 = np.cross(s2, s0)
    e2 = np.cross(s2, e1)
    e1 /= norm(e1)
    e2 /= norm(e2)
    e3 = s2 / norm(s2)
    R = np.array([e1, e2, e3], dtype=np.float64)
    return R
