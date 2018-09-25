#
# See https://github.com/dials/dials/wiki/pytest for documentation on how to
# write and run pytest tests, and an overview of the available features.
#

from __future__ import absolute_import, division, print_function

import os

import py
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
def xia2_regression():
  '''Return the absolute path to the xia2_regression module as a string.
     Skip the test if xia2_regression is not installed.'''
  try:
    import xia2_regression as xr
  except ImportError:
    pytest.skip("xia2_regression required for this test")
  return os.path.dirname(xr.__file__)

@pytest.fixture(scope="session")
def xia2_regression_build():
  '''Return the absolute path to the xia2_regression directory within the build
     path as a string. Skip the test if xia2_regression is not installed.'''
  try:
    import xia2_regression as xr
  except ImportError:
    pytest.skip("xia2_regression required for this test")
  try:
    dls_common_copy = '/dls/science/groups/scisoft/DIALS/repositories/current/xia2_regression_data'
    if os.path.exists(dls_common_copy):
      return dls_common_copy
  except Exception:
    pass
  try:
    x2rpath = os.path.join(os.environ.get('LIBTBX_BUILD'), 'xia2_regression')
  except AttributeError:
    x2rpath = ''
  if not os.path.exists(x2rpath):
    pytest.skip("xia2_regression required for this test")
  if 'test_data' not in os.listdir(x2rpath):
    pytest.skip("xia2_regression files need to be downloaded for this test. Run xia2_regression.fetch_test_data")
  return x2rpath

@pytest.fixture(scope="session")
def regression_data():
  '''Return the location of a regression data set as py.path object.
     Skip the test if the data are not present.
  '''
  dls_dir = '/dls/science/groups/scisoft/DIALS/regression_data'
  if os.getenv('REGRESSIONDATA'):
    target_dir = os.getenv('REGRESSIONDATA')
  elif os.path.exists(os.path.join(dls_dir, 'filelist.json')):
    target_dir = dls_dir
  elif os.getenv('LIBTBX_BUILD'):
    target_dir = os.path.join(os.getenv('LIBTBX_BUILD'), 'regression_data')
  else:
    pytest.skip('Can not determine regression data location. Use environment variable REGRESSIONDATA')

  from xia2.Test.fetch_test_data import download_lock, fetch_test_data
  class DataFetcher():
    _cache = {}
    def __call__(self, test_data):
      if test_data not in self._cache:
        with download_lock(target_dir):
          self._cache[test_data] = fetch_test_data(target_dir, pre_scan=True, file_group=test_data, read_only=True)
          if self._cache[test_data]:
            self._cache[test_data] = str(self._cache[test_data]) # https://github.com/cctbx/cctbx_project/issues/234
      if not self._cache[test_data]:
        pytest.skip('Regression data is required to run this test. Run xia2.fetch_test_data')
      return py.path.local(self._cache[test_data])
    def __repr__(self):
      return "<R/O DataFetcher: %s>" % target_dir

  return DataFetcher()

@pytest.fixture
def run_in_tmpdir(tmpdir):
  '''Shortcut to create a temporary directory and then run the test inside
     this directory.'''
  tmpdir.chdir()
  return tmpdir
