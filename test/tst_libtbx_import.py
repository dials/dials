from __future__ import division
import libtbx.load_env
from libtbx import easy_run
import os

def run():
  cmd = os.path.join(libtbx.env.under_build('bin'), 'libtbx.import_all_ext')
  easy_run.fully_buffered(cmd).raise_if_errors()
  print "OK"

if __name__ == '__main__':
  run()
