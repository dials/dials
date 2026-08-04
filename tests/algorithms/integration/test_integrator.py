from __future__ import annotations

import pickle
from unittest import mock

from dxtbx.model import Beam, Crystal, Experiment, ExperimentList, Goniometer, Scan

from dials.algorithms.integration import integrator


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


def test_frame_sliced_integrator_executor_is_picklable():
    # A real experiment list is required, as the executor calculates the
    # frame orientations on construction
    experiments = ExperimentList(
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
            )
        ]
    )
    executor = integrator.FrameSlicedIntegratorExecutor(
        experiments=experiments,
    )
    assert len(executor.frame_orientations) == len(experiments)
    assert len(executor.frame_orientations[0].UB) == 5

    pickled = pickle.dumps(executor)
    unpickled = pickle.loads(pickled)
    assert isinstance(unpickled, integrator.FrameSlicedIntegratorExecutor)
    assert len(unpickled.frame_orientations) == len(experiments)
