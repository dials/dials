#
# See https://github.com/dials/dials/wiki/pytest for documentation on how to
# write and run pytest tests, and an overview of the available features.
#


from __future__ import annotations

import multiprocessing
import os
import sys
import warnings
from pathlib import Path

import pytest

if sys.version_info[:2] == (3, 7) and sys.platform == "darwin":
    multiprocessing.set_start_method("forkserver")

collect_ignore = []


def _build_filterwarnings_string() -> str:
    """
    Given the set of active warnings, build a PYTHONWARNINGS string.

    This lets us set the PYTHONWARNINGS environment variable so that
    warning filters are passed down to subprocesses.
    """
    filter_parts = []
    for action, regex, category, modregex, line in warnings.filters:
        if action != "ignore":
            continue
        if category.__module__ != "builtins":
            continue
        this_action = [
            action,
            regex.pattern if regex else "",
            category.__name__,
        ]
        if modregex is not None:
            this_action.append(modregex.pattern)
            if line:
                this_action.append(line)
        filter_parts.append(":".join(str(x) for x in this_action))

    return ",".join(filter_parts)


def pytest_configure(config):
    if not config.pluginmanager.hasplugin("dials_data"):

        @pytest.fixture(scope="session")
        def dials_data():
            pytest.skip("This test requires the dials_data package to be installed")

        globals()["dials_data"] = dials_data

    # Ensure that subprocesses get the warnings filters
    os.environ["PYTHONWARNINGS"] = _build_filterwarnings_string()


@pytest.fixture(scope="session")
def dials_regression() -> Path:
    """Return the absolute path to the dials_regression module as a string.
    Skip the test if dials_regression is not installed."""

    if "DIALS_REGRESSION" in os.environ:
        return os.environ["DIALS_REGRESSION"].rstrip("/")

    try:
        import dials_regression as dr

        return Path(dr.__file__).parent
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
            return Path(reference_copy)
    except ImportError:
        pass  # Cannot tell whether in DLS network or not
    pytest.skip("dials_regression required for this test")


@pytest.fixture
def run_in_tmp_path(tmp_path) -> Path:
    """
    A fixture to change the working directory for the test to a temporary directory.

    The original working directory is restored upon teardown of the fixture.

    Args:
        tmp_path: Pytest tmp_path fixture, see
                  https://docs.pytest.org/en/latest/how-to/tmp_path.html

    Yields:
        The path to the temporary working directory defined by tmp_path.
    """
    cwd = Path.cwd()
    os.chdir(tmp_path)
    yield tmp_path
    os.chdir(cwd)
