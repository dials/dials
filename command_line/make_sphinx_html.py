# LIBTBX_SET_DISPATCHER_NAME dev.dials.make_sphinx_html

from __future__ import division
from dials.util.procrunner import run_process
import libtbx.load_env
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

if __name__ == "__main__":
  cctbx_base = libtbx.env.find_in_repositories("cctbx_project")
  dials_dir = libtbx.env.find_in_repositories("dials")
  dials_github_io = libtbx.env.find_in_repositories("dials.github.io")
  assert (dials_github_io is not None)
  assert (cctbx_base is not None)
  base_dir = os.path.dirname(cctbx_base)
  dest_dir = dials_github_io
  os.chdir(os.path.join(dials_dir, "doc", "sphinx"))
  result = run_process(["make", "clean"])
  assert result['exitcode'] == 0, \
      'make clean failed with exit code %d' % result['exitcode']
  result = run_process(["make", "html"])
  assert result['exitcode'] == 0, \
      'make html failed with exit code %d' % result['exitcode']
  print "Copying HTML pages to", dest_dir
  recursive_overwrite("build/html", dest_dir)
