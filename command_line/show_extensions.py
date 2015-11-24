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

# LIBTBX_SET_DISPATCHER_NAME dev.dials.show_extensions

class Script(object):
  ''' The class to encapsulate the script. '''

  def __init__(self):
    ''' Initialise the script. '''
    from dials.util.options import OptionParser
    from libtbx.phil import parse
    import libtbx.load_env

    # Create the phil parameters
    phil_scope = parse('''
      interfaces=False
        .type = bool
        .help = "Only show information about the interfaces"
    ''')

    # Create the option parser
    usage = "usage: %s [options] /path/to/image/files" \
      % libtbx.env.dispatcher_name
    self.parser = OptionParser(usage=usage, phil=phil_scope)

  def run(self):
    ''' Run the script. '''
    import dials.extensions # import dependency
    from dials.interfaces import ProfileModelIface
    from dials.interfaces import BackgroundIface
    from dials.interfaces import CentroidIface
    from dials.interfaces import SpotFinderThresholdIface

    # Parse the command line arguments
    params, options = self.parser.parse_args()

    # Create the list of interfaces
    interfaces = [
      ProfileModelIface,
      BackgroundIface,
      CentroidIface,
      SpotFinderThresholdIface
    ]

    # Loop through all the interfaces
    for iface in interfaces:
      print '-' * 80
      print 'Interface: %s' % iface.__name__

      # Either just show information about interfaces or show some about
      # extensions depending on user input
      if params.interfaces:

        # Print info about interface
        if options.verbose > 0:
          print ' name = %s' % iface.name
          if options.verbose > 1:
            level = options.verbose - 2
            scope = iface.phil_scope()
            phil = scope.as_str(print_width=80-3, attributes_level=level)
            phil = '\n'.join((' ' * 2) + l for l in phil.split('\n'))
            if phil.strip() != '':
              print ' phil:\n%s' % phil

      else:

        # Loop through all the extensions
        for ext in iface.extensions():
          print ' Extension: %s' % ext.__name__
          if options.verbose > 0:
            print '  name = %s' % ext.name
            if options.verbose > 1:
              level = options.verbose - 2
              scope = ext.phil_scope()
              phil = scope.as_str(print_width=80-3, attributes_level=level)
              phil = '\n'.join((' ' * 3) + l for l in phil.split('\n'))
              if phil.strip() != '':
                print '  phil:\n%s' % phil


if __name__ == '__main__':
  from dials.util import halraiser
  try:
    script = Script()
    script.run()
  except Exception as e:
    halraiser(e)
