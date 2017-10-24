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
  try:
    return libtbx.env.dist_path('dials_regression')
  except KeyError:
    pytest.skip("dials_regression required for this test")

@pytest.fixture
def xia2_regression():
  try:
    return libtbx.env.dist_path('xia2_regression')
  except KeyError:
    pytest.skip("xia2_regression required for this test")

@pytest.fixture
def xia2_regression_build():
  x2rpath = libtbx.env.under_build('xia2_regression')
  if not os.path.exists(x2rpath):
    pytest.skip("xia2_regression required for this test")
  return x2rpath

from libtbx.test_utils.pytest import libtbx_collector
pytest_collect_file = libtbx_collector()
