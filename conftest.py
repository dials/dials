#
# See https://github.com/dials/dials/wiki/pytest for documentation on how to
# write and run pytest tests, and an overview of the available features.
#

from __future__ import absolute_import, division, print_function

import os

import pytest

try:
    import dials_data as pkg_dials_data

    dials_data = pkg_dials_data.dials_data
except ImportError:
    pkg_dials_data = None

    @pytest.fixture
    def dials_data():
        pytest.skip("Test requires python package dials_data")


def pytest_addoption(parser):
    """Add a '--runslow' option to py.test."""
    if pkg_dials_data:
        pkg_dials_data.pytest_addoption(parser)

    try:
        parser.addoption(
            "--regression",
            action="store_true",
            default=False,
            help="run regression tests",
        )
    except ValueError:
        pass  # Thrown in case the command line option is already defined
    parser.addoption(
        "--runslow", action="store_true", default=False, help="run slow tests"
    )


def pytest_collection_modifyitems(config, items):
  '''Tests marked as slow will not be run unless slow tests are enabled with
     the '--runslow' parameter or the test is selected specifically. The
     latter allows running slow tests via the libtbx compatibility layer.'''
  if not config.getoption("--runslow") and len(items) > 1:
    skip_slow = pytest.mark.skip(reason="need --runslow option to run")
    for item in items:
      if "slow" in item.keywords:
        item.add_marker(skip_slow)

@pytest.fixture(scope="session")
def dials_regression():
  '''Return the absolute path to the dials_regression module as a string.
     Skip the test if dials_regression is not installed.'''
  try:
    import dials_regression as dr
  except ImportError:
    pytest.skip("dials_regression required for this test")
  return os.path.dirname(dr.__file__)


@pytest.fixture
def run_in_tmpdir(tmpdir):
  '''Shortcut to create a temporary directory and then run the test inside
     this directory.'''
  tmpdir.chdir()
  return tmpdir
