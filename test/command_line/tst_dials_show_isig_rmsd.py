from __future__ import division
import os
from libtbx import easy_run
from libtbx.test_utils import open_tmp_directory

import libtbx.load_env
have_dials_regression = libtbx.env.has_module("dials_regression")
if have_dials_regression:
  dials_regression = libtbx.env.find_in_repositories(
    relative_path="dials_regression",
    test=os.path.isdir)
else:
  dials_regression = None

def exercise1():
  if not dials_regression:
    return
  integrated_pickle = os.path.join(
    dials_regression, 'centroid_test_data', 'integrated.pickle')
  cwd = os.path.abspath(os.curdir)
  tmp_dir = open_tmp_directory()
  os.chdir(tmp_dir)
  cmd = " ".join(["dials.show_isig_rmsd", integrated_pickle])
  result = easy_run.fully_buffered(cmd).raise_if_errors()
  if result:
    print 'OK'
  os.chdir(cwd)

def run(args):
  if not have_dials_regression:
    print "Skipping tst_dials_show_isig_rmsd.py: dials_regression not available"
    return

  exercise1()

if __name__ == '__main__':
  import sys
  from libtbx.utils import show_times_at_exit
  show_times_at_exit()
  run(sys.argv[1:])
