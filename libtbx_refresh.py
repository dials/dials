from __future__ import division

def run():
  from dials.util.config import CompletionGenerator
  gen = CompletionGenerator()
  # generate init.sh and SConscript file in build/dials/autocomplete
  gen.generate()
# # add init.sh to setpaths.sh
# gen.install() 

run()

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

try:
  from dials.framework import env
  env.cache.wipe()
except Exception:
  pass
