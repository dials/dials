from __future__ import division

from math import exp, sqrt

from scitbx import matrix
from scitbx.linalg import eigensystem, l_l_transpose_cholesky_decomposition_in_place

from dials.algorithms.profile_model.potato import (
    BBoxCalculatorAngular,
    BBoxCalculatorSimple,
    MaskCalculatorAngular,
    MaskCalculatorSimple,
    PredictorAngular,
    PredictorSimple,
)
from dials.algorithms.profile_model.potato.parameterisation import (
    Angular2MosaicityParameterisation,
    Angular4MosaicityParameterisation,
    Simple1MosaicityParameterisation,
    Simple6MosaicityParameterisation,
)
from dials.array_family import flex


class ProfileModelBase(object):
    """
    Class to store profile model

    """

    def __init__(self, params):
        """
        Initialise the class

        """
        self.params = params

    def sigma(self):
        """
        Get the sigma

        """
        return self.parameterisation().sigma()

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
        eigen_decomposition = eigensystem.real_symmetric(
            state.mosaicity_covariance_matrix.as_flex_double_matrix()
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
            "__id__": self.__class__.__name__,
            "parameters": params,
            "sigma": sigma.as_numpy_array().tolist(),
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

    def predict_reflections(self, experiments, miller_indices, probability=0.9973):
        """
        Predict the reflections

        """
        predictor = PredictorSimple(experiments[0], self.sigma(), probability)
        return predictor.predict(miller_indices)

    def compute_bbox(self, experiments, reflections, probability=0.9973):
        """
        Compute the bounding box

        """
        calculator = BBoxCalculatorSimple(experiments[0], self.sigma(), probability, 4)
        calculator.compute(reflections)

    def compute_mask(self, experiments, reflections, probability=0.9973):
        """
        Compute the mask

        """
        calculator = MaskCalculatorSimple(experiments[0], self.sigma(), probability)
        calculator.compute(reflections)

    def sigma_for_reflection(self, s0, r):
        """
        Get sigma for a reflections

        """
        return self.sigma()

    def compute_partiality(self, experiments, reflections):
        """
        Compute the partiality

        """
        s0 = matrix.col(experiments[0].beam.get_s0())
        # num = reflections.get_flags(reflections.flags.indexed).count(True)
        num = reflections.size()

        # Compute the marginal variance for the 000 reflection
        S00 = experiments[0].crystal.mosaicity.sigma()[8]

        partiality = flex.double(len(reflections))
        partiality_variance = flex.double(len(reflections))
        for k in range(len(reflections)):
            s2 = matrix.col(reflections[k]["s2"])
            sigma = experiments[0].crystal.mosaicity.sigma()
            R = compute_change_of_basis_operation(s0, s2)
            S = R * (sigma) * R.transpose()
            mu = R * s2
            assert abs(1 - mu.normalize().dot(matrix.col((0, 0, 1)))) < 1e-7
            S22 = S[8]
            mu2 = mu[2]
            eps = s0.length() - mu2
            var_eps = S22 / num  # FIXME Approximation
            partiality[k] = exp(-0.5 * eps * (1 / S22) * eps) * sqrt(S00 / S22)
            partiality_variance[k] = (
                var_eps * (eps ** 2 / (S00 * S22)) * exp(eps ** 2 / S22)
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

    name = "Simple1ProfileModel"

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
        return Class.from_params(flex.double((sigma_d,)))

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

    name = "Simple6ProfileModel"

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
        return Class.from_params(flex.double((sigma_d, 0, sigma_d, 0, 0, sigma_d)))

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
        return Q.transpose() * self.sigma() * Q

    def predict_reflections(self, experiments, miller_indices, probability=0.9973):
        """
        Predict the reflections

        """
        predictor = PredictorAngular(experiments[0], self.sigma(), probability)
        return predictor.predict(miller_indices)

    def compute_bbox(self, experiments, reflections, probability=0.9973):
        """
        Compute the bounding box

        """
        calculator = BBoxCalculatorAngular(experiments[0], self.sigma(), probability, 4)
        calculator.compute(reflections)

    def compute_mask(self, experiments, reflections, probability=0.9973):
        """
        Compute the mask

        """
        calculator = MaskCalculatorAngular(experiments[0], self.sigma(), probability)
        calculator.compute(reflections)

    def compute_partiality(self, experiments, reflections):
        """
        Compute the partiality

        """
        s0 = matrix.col(experiments[0].beam.get_s0())
        num = reflections.get_flags(reflections.flags.indexed).count(True)
        num = reflections.size()
        partiality = flex.double(len(reflections))
        partiality_variance = flex.double(len(reflections))
        for k in range(len(reflections)):
            s2 = matrix.col(reflections[k]["s2"])
            r = s2 - s0
            sigma = experiments[0].crystal.mosaicity.sigma()
            R = compute_change_of_basis_operation(s0, s2)
            Q = compute_change_of_basis_operation(s0, r)
            S = R * (Q.transpose() * sigma * Q) * R.transpose()
            mu = R * s2
            assert abs(1 - mu.normalize().dot(matrix.col((0, 0, 1)))) < 1e-7
            S22 = S[8]
            mu2 = mu[2]
            eps = s0.length() - mu2
            var_eps = S22 / num  # FIXME Approximation
            S00 = S22  # FIXME
            partiality[k] = exp(-0.5 * eps * (1 / S22) * eps) * sqrt(S00 / S22)
            partiality_variance[k] = (
                var_eps * (eps ** 2 / (S00 * S22)) * exp(eps ** 2 / S22)
            )

        reflections["partiality"] = partiality
        reflections["partiality.inv.variance"] = partiality_variance

    @classmethod
    def from_params(Class, params):
        """
        Create the class from some parameters

        """
        return Class(params)


class Angular2ProfileModel(AngularProfileModelBase):
    """
    Class to store profile model

    """

    name = "Angular2ProfileModel"

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
        return Class.from_params(flex.double((sigma_d, sigma_d)))

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

    name = "Angular4ProfileModel"

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
        return Class.from_params(flex.double((sigma_d, 0, sigma_d, sigma_d)))

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
        if model == "simple1":
            return Simple1ProfileModel.from_sigma_d(sigma_d)
        elif model == "simple6":
            return Simple6ProfileModel.from_sigma_d(sigma_d)
        elif model == "angular2":
            return Angular2ProfileModel.from_sigma_d(sigma_d)
        elif model == "angular4":
            return Angular4ProfileModel.from_sigma_d(sigma_d)

        raise RuntimeError(f"Unknown profile model: {model}")


def compute_change_of_basis_operation(s0, s2):
    """
    Compute the change of basis operation that puts the s2 vector along the z axis

    """
    e1 = s2.cross(s0).normalize()
    e2 = s2.cross(e1).normalize()
    e3 = s2.normalize()
    R = matrix.sqr(e1.elems + e2.elems + e3.elems)
    return R
