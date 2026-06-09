from __future__ import annotations

import pickle
from unittest import mock

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
