# LIBTBX_SET_DISPATCHER_NAME dev.dials.make_sphinx_html

from __future__ import division
from libtbx import easy_run
import libtbx.load_env
import os.path as op
import shutil
import os

def recursive_overwrite(src, dest, ignore=None):
  if os.path.isdir(src):
    if not os.path.isdir(dest):
      os.makedirs(dest)
    files = os.listdir(src)
    if ignore is not None:
      ignored = ignore(src, files)
    else:
      ignored = set()
    for f in files:
      if f not in ignored:
        recursive_overwrite(os.path.join(src, f),
                            os.path.join(dest, f),
                            ignore)
  else:
    shutil.copyfile(src, dest)

if (__name__ == "__main__") :
  cctbx_base = libtbx.env.find_in_repositories("cctbx_project")
  dials_dir = libtbx.env.find_in_repositories("dials")
  dials_htdocs = libtbx.env.find_in_repositories("dials_htdocs")
  if dials_htdocs is None:
    dials_htdocs = libtbx.env.find_in_repositories("htdocs")
  assert (dials_htdocs is not None)
  assert (cctbx_base is not None)
  base_dir = op.dirname(cctbx_base)
  dest_dir = dials_htdocs
  os.chdir(op.join(dials_dir, "doc", "sphinx"))
  easy_run.call("make clean")
  easy_run.call("make html")
  print "Moving HTML pages to", dest_dir
  recursive_overwrite("build/html", dest_dir)
