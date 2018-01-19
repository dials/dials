#
# See https://github.com/dials/dials/wiki/pytest for documentation on how to
# write and run pytest tests, and an overview of the available features.
#

from __future__ import absolute_import, division, print_function

import libtbx.load_env
import os
import pytest

@pytest.fixture
def dials_regression():
  '''Return the absolute path to the dials_regression module as a string.
     Skip the test if dials_regression is not installed.'''
  try:
    import dials_regression as dr
  except ImportError:
    pytest.skip("dials_regression required for this test")
  return os.path.dirname(dr.__file__)

@pytest.fixture
def xia2_regression():
  '''Return the absolute path to the xia2_regression module as a string.
     Skip the test if dials_regression is not installed.'''
  try:
    import xia2_regression as xr
  except ImportError:
    pytest.skip("xia2_regression required for this test")
  return os.path.dirname(xr.__file__)

@pytest.fixture
def xia2_regression_build():
  '''Return the absolute path to the xia2_regression directory within the build
     path as a string. Skip the test if xia2_regression is not installed.'''
  x2rpath = libtbx.env.under_build('xia2_regression')
  if not os.path.exists(x2rpath):
    pytest.skip("xia2_regression required for this test")
  return x2rpath

from libtbx.test_utils.pytest import libtbx_collector
pytest_collect_file = libtbx_collector()
