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
  filenames = glob("extensions/*.pyc")
  if len(filenames) > 0:
    print "Cleaning up 'dials/extensions':"
    for filename in glob("extensions/*.pyc"):
      print " Deleting %s" % filename
      os.remove(filename)
except Exception, e:
  pass
