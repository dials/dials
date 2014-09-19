from __future__ import division

def run():
  from dials.util.config import CompletionGenerator
  gen = CompletionGenerator()
  gen.generate()

try:
  run()
except Exception, e:
  pass

try:
  from glob import glob
  import os
  from os.path import join
  import libtbx.load_env
  dials_path = libtbx.env.dist_path('dials')
  filenames = glob(join(dials_path, "extensions", "*.pyc"))
  if len(filenames) > 0:
    print "Cleaning up 'dials/extensions':"
    for filename in filenames:
      print " Deleting %s" % filename
      os.remove(filename)
except Exception, e:
  pass
