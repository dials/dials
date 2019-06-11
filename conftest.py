#
# See https://github.com/dials/dials/wiki/pytest for documentation on how to
# write and run pytest tests, and an overview of the available features.
#

from __future__ import absolute_import, division, print_function

import os

import pytest
import six


def pytest_addoption(parser):
    """Add a '--runslow' option to py.test."""
    parser.addoption(
        "--runslow", action="store_true", default=False, help="run slow tests"
    )


def pytest_collection_modifyitems(config, items):
    """Tests marked as slow will not be run unless slow tests are enabled with
    the '--runslow' parameter or the test is selected specifically. The
    latter allows running slow tests via the libtbx compatibility layer."""
    if not config.getoption("--runslow") and len(items) > 1:
        skip_slow = pytest.mark.skip(reason="need --runslow option to run")
        for item in items:
            if "slow" in item.keywords:
                item.add_marker(skip_slow)


def pytest_configure(config):
    if six.PY3:
        import dxtbx.tests.python3_test_filter as ptf

        exp = ptf.Python3TestFailureExpectationPlugin(config)
        config.pluginmanager.register(exp)


@pytest.fixture(scope="session")
def dials_regression():
    """Return the absolute path to the dials_regression module as a string.
    Skip the test if dials_regression is not installed."""
    try:
        import dials_regression as dr

        return os.path.dirname(dr.__file__)
    except ImportError:
        pass  # dials_regression not configured
    try:
        import socket

        reference_copy = "/dls/science/groups/scisoft/DIALS/repositories/git-reference/dials_regression"
        if (
            os.name == "posix"
            and "diamond.ac.uk" in socket.gethostname()
            and os.path.exists(reference_copy)
        ):
            return reference_copy
    except ImportError:
        pass  # can not tell whether in DLS network or not
    pytest.skip("dials_regression required for this test")


@pytest.fixture
def run_in_tmpdir(tmpdir):
    """Shortcut to create a temporary directory and then run the test inside
    this directory."""
    tmpdir.chdir()
    return tmpdir
