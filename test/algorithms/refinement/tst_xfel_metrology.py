#!/usr/bin/env cctbx.python

#
#  Copyright (C) (2016) STFC Rutherford Appleton Laboratory, UK.
#
#  Author: David Waterman.
#
#  This code is distributed under the BSD license, a copy of which is
#  included in the root directory of this package.
#

"""
A basic test of joint refinement of the CS-PAD detector at hierarchy level 2
with 300 crystals.

"""

# python imports
from __future__ import division
import os
import shutil
import libtbx.load_env # required for libtbx.env.find_in_repositories
from libtbx import easy_run
from libtbx.test_utils import open_tmp_directory
from dials.array_family import flex
import cPickle as pickle

class Test(object):

  def __init__(self):

    dials_regression = libtbx.env.find_in_repositories(
      relative_path="dials_regression",
      test=os.path.isdir)

    self._data_dir = os.path.join(dials_regression, "refinement_test_data",
                            "xfel_metrology")

    # set up a temporary directory
    self._cwd = os.path.abspath(os.curdir)
    self._tmp_dir = open_tmp_directory()
    os.chdir(self._tmp_dir)

  def _finish(self):
    """Change back to the original directory and delete the temp dir"""

    os.chdir(self._cwd)
    shutil.rmtree(self._tmp_dir)
    return

  def _refine(self):
    """Do refinement and load the history"""

    experiments_path = os.path.join(self._data_dir, 'benchmark_level2d.json')
    reflections_path = os.path.join(self._data_dir, 'benchmark_level2d.pickle')
    phil_path = os.path.join(self._data_dir, 'refine.phil')
    cmd="dials.refine {0} {1} {2} history=history.pickle".format(
      experiments_path, reflections_path, phil_path)
    result = easy_run.fully_buffered(command=cmd).raise_if_errors()

    # there are plenty of things we could do with the refinement history, but
    # here just check that final RMSDs are low enough
    with open('history.pickle', 'r') as f:
      self._history = pickle.load(f)
    final_rmsd = self._history['rmsd'][-1]
    assert final_rmsd[0] < 0.0354
    assert final_rmsd[1] < 0.0406
    assert final_rmsd[2] < 0.0018

    return

  def run(self):
    """Run the test"""

    try:
      self._refine()
    finally:
      self._finish()

    return

def run():
  if not libtbx.env.has_module("dials_regression"):
    print "Skipping tests in " + __file__ + " as dials_regression not present"
    return
  try:
    from scitbx.examples.bevington import non_linear_ls_eigen_wrapper
  except ImportError:
    print "Skipping tests in " + __file__ + " as SparseLevMar engine not available"

  tst = Test()
  tst.run()

if __name__ == '__main__':
  from libtbx.utils import show_times_at_exit
  show_times_at_exit()
  run()
  print "OK"
