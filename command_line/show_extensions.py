#!/usr/bin/env python
#
# show_extensions.py
#
#  Copyright (C) 2013 Diamond Light Source
#
#  Author: James Parkhurst
#
#  This code is distributed under the BSD license, a copy of which is
#  included in the root directory of this package.
from __future__ import division

if __name__ == '__main__':
  from optparse import OptionParser
  from dials.framework.registry import Registry

  usage = "usage: %prog [options] /path/to/image/files"
  parser = OptionParser(usage)

  # Print verbose output
  parser.add_option(
    "-v", "--verbose",
    dest = "verbose",
    action = "count", default = 0,
    help = "Set the verbosity level (-vv gives a verbosity level of 2)")

  # Parse the command line arguments
  options, args = parser.parse_args()

  # Get the registry
  registry = Registry()

  # Loop through all the interfaces
  for iface in registry.interfaces():
    print '-' * 80
    print 'Interface: %s' % iface.__name__

    # Loop through all the extensions
    for ext in iface.extensions():
      print ' Extension: %s' % ext.__name__
      if options.verbose > 0:
        print '  name = %s' % ext.name
        if options.verbose > 1:
          level = options.verbose - 2
          phil = ext.phil_scope().as_str(print_width=80-4, attributes_level=level)
          phil = '\n'.join((' ' * 4) + l for l in phil.split('\n'))
          print '  phil:\n%s' % phil
