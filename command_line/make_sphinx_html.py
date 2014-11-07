# LIBTBX_SET_DISPATCHER_NAME dev.dials.make_sphinx_html

from __future__ import division
from libtbx import easy_run
import libtbx.load_env
import os.path as op
import shutil
import os

if (__name__ == "__main__") :
  cctbx_base = libtbx.env.find_in_repositories("cctbx_project")
  dials_dir = libtbx.env.find_in_repositories("dials")
  dials_htdocs = libtbx.env.find_in_repositories("dials_htdocs")
  if dials_htdocs is None:
    dials_htdocs = libtbx.env.find_in_repositories("htdocs")
  assert (dials_htdocs is not None)
  assert (cctbx_base is not None)
  base_dir = op.dirname(cctbx_base)
  dest_dir = op.join(dials_htdocs, "doc")
  if op.exists(dest_dir) :
    shutil.rmtree(dest_dir)
  os.chdir(op.join(dials_dir, "doc", "sphinx"))
  easy_run.call("make html")
  shutil.move("build/html", dest_dir)
