#
# See https://github.com/dials/dials/wiki/pytest for documentation on how to
# write and run pytest tests, and an overview of the available features.
#


from __future__ import absolute_import, division, print_function

import os
import warnings

import pytest
import six

# https://stackoverflow.com/a/40846742
warnings.filterwarnings("ignore", message="numpy.dtype size changed")
warnings.filterwarnings("ignore", message="numpy.ufunc size changed")

collect_ignore = []
if six.PY2:
    _base = os.path.dirname(__file__)
    with open(os.path.join(_base, ".travis", "python2-supported-files"), "r") as fh:
        allowed_testfiles = {tuple(f.strip()[2:].split("/")) for f in fh}
    for root, dirs, files in os.walk(_base):
        relroot = os.path.relpath(root, _base).split(os.path.sep)
        if relroot == ["."]:
            relroot = []
        for f in files:
            if f.endswith(".py"):
                filetuple = tuple(relroot + [f])
                if filetuple not in allowed_testfiles:
                    collect_ignore.append(os.path.join(*filetuple))
    warnings.warn(
        "%d test files were excluded as they can only be interpreted with Python 3"
        % len(collect_ignore),
        UserWarning,
    )


@pytest.fixture(scope="session")
def python3():
    if six.PY2:
        pytest.skip("Test requires a Python 3 installation")


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
    if not config.pluginmanager.hasplugin("dials_data"):

        @pytest.fixture(scope="session")
        def dials_data():
            pytest.skip("This test requires the dials_data package to be installed")

        globals()["dials_data"] = dials_data


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
            and socket.gethostname().endswith(".diamond.ac.uk")
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
