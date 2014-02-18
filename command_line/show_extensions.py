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
  from dials.framework.registry import Registry

  # Get the registry
  registry = Registry()

  # Loop through all the interfaces
  for iface in registry.interfaces():
    print '-' * 80
    print 'Interface: %s' % iface.__name__

    # Loop through all the extensions
    for ext in iface.extensions():
      print ' Extension: %s' % ext.__name__
