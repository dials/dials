from __future__ import division

def version():
  from dials.util.version import dials_version
  import dials
  import os

  print dials_version()
  print "Installed in: %s" % os.path.split(dials.__file__)[0]

version()
