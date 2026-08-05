from __future__ import annotations

import pickle
from unittest import mock

import pytest

from dxtbx.model import Beam, Crystal, Experiment, ExperimentList, Goniometer, Scan

from dials.algorithms.integration import integrator
from dials.algorithms.profile_model.gaussian_rs import Model as GaussianRSProfileModel
from dials.array_family import flex


def test_profile_modeller_executor_is_picklable():
    executor = integrator.ProfileModellerExecutor(
        experiments=mock.ANY,
        profile_fitter=mock.ANY,
    )
    pickled = pickle.dumps(executor)
    unpickled = pickle.loads(pickled)
    assert isinstance(unpickled, integrator.ProfileModellerExecutor)


def test_profile_validator_executor_is_picklable():
    executor = integrator.ProfileValidatorExecutor(
        experiments=mock.ANY,
        profile_fitter=mock.ANY,
    )
    pickled = pickle.dumps(executor)
    unpickled = pickle.loads(pickled)
    assert isinstance(unpickled, integrator.ProfileValidatorExecutor)


def test_integrator_executor_is_picklable():
    executor = integrator.IntegratorExecutor(
        experiments=mock.ANY,
    )
    pickled = pickle.dumps(executor)
    unpickled = pickle.loads(pickled)
    assert isinstance(unpickled, integrator.IntegratorExecutor)


def frame_sliced_integrator_executor_experiments(sigma_m=0.1):
    """A single experiment with a five image scan, which is all the frame sliced
    executor needs to construct itself"""
    return ExperimentList(
        [
            Experiment(
                beam=Beam((0, 0, 1), 1.0),
                goniometer=Goniometer((1, 0, 0)),
                scan=Scan((1, 5), (0.0, 0.5)),
                crystal=Crystal(
                    real_space_a=(10, 0, 0),
                    real_space_b=(0, 10, 0),
                    real_space_c=(0, 0, 10),
                    space_group_symbol="P 1",
                ),
                profile=GaussianRSProfileModel(
                    params=None, n_sigma=3, sigma_b=0.02, sigma_m=sigma_m, deg=False
                ),
            )
        ]
    )


def test_frame_sliced_integrator_executor_is_picklable():
    # A real experiment list is required, as the executor calculates the frame
    # orientations on construction and takes the mosaic spread from the profile
    # model
    experiments = frame_sliced_integrator_executor_experiments()
    executor = integrator.FrameSlicedIntegratorExecutor(
        experiments=experiments,
    )
    assert len(executor.frame_orientations) == len(experiments)
    assert len(executor.frame_orientations[0].UB) == 5

    pickled = pickle.dumps(executor)
    unpickled = pickle.loads(pickled)
    assert isinstance(unpickled, integrator.FrameSlicedIntegratorExecutor)
    assert len(unpickled.frame_orientations) == len(experiments)


def test_frame_sliced_integrator_executor_rocking_curve_width():
    """The width of the rocking curve in the rotation angle is sigma_m / |zeta|,
    with the sign of zeta making no difference"""
    executor = integrator.FrameSlicedIntegratorExecutor(
        experiments=frame_sliced_integrator_executor_experiments(sigma_m=0.1),
    )

    zeta = flex.double([1.0, 0.5, -0.5])
    z_cal = flex.double([0.5, 1.5, 2.5])
    assert list(executor._sigma_phi(zeta, z_cal, 0)) == pytest.approx([0.1, 0.2, 0.2])


def test_frame_sliced_integrator_executor_scan_varying_rocking_curve_width():
    """A scan-varying mosaic spread is looked up at the predicted centroid frame
    of the reflection, as PartialityCalculator3D does"""
    executor = integrator.FrameSlicedIntegratorExecutor(
        experiments=frame_sliced_integrator_executor_experiments(
            sigma_m=flex.double([0.1, 0.2, 0.3, 0.4, 0.5])
        ),
    )

    # The five images of the scan cover array indices 0 to 5
    zeta = flex.double(7, 1.0)
    z_cal = flex.double([0.0, 0.5, 1.5, 2.9, 4.5, -0.5, 5.5])
    # The first and last are off the ends of the scan, so take the mosaic spread
    # of the first and last images
    assert list(executor._sigma_phi(zeta, z_cal, 0)) == pytest.approx(
        [0.1, 0.1, 0.2, 0.3, 0.5, 0.1, 0.5]
    )
