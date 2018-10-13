#
# See https://github.com/dials/dials/wiki/pytest for documentation on how to
# write and run pytest tests, and an overview of the available features.
#

from __future__ import absolute_import, division, print_function

import os

import pytest

def pytest_addoption(parser):
  '''Add a '--runslow' option to py.test.'''
  parser.addoption("--runslow", action="store_true", default=False,
                   help="run slow tests")

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

@pytest.fixture(scope="session")
def regression_data():
  '''Return the location of a regression data set as py.path object.
     Skip the test if the data are not present.
  '''
  import dials.util.regression_data
  df = dials.util.regression_data.DataFetcher(read_only=True)
  def skip_test_if_lookup_failed(result):
    if not result:
      pytest.skip('Regression data is required to run this test. Run dials.fetch_test_data')
    return result
  setattr(df, 'result_filter', skip_test_if_lookup_failed)
  return df

@pytest.fixture
def run_in_tmpdir(tmpdir):
  '''Shortcut to create a temporary directory and then run the test inside
     this directory.'''
  tmpdir.chdir()
  return tmpdir
