#!/usr/bin/env python
#
# h5dump.py
#
#  Copyright (C) 2013 Diamond Light Source
#
#  Author: James Parkhurst
#
#  This code is distributed under the BSD license, a copy of which is
#  included in the root directory of this package.
#
from __future__ import division

if __name__ == '__main__':
  from dials.array_family import flex
  from optparse import OptionParser
  from dials.util.command_line import Command

  # Specify the command line options
  usage  = "usage: %prog [options] " \
           "/path/to/reflection/file.pkl " \
           "/path/to/reflection/file.h5"

  # Create an option parser
  parser = OptionParser(usage)

  # Parse the arguments
  options, args = parser.parse_args()

  if len(args) < 2:
    parser.print_help()

  else:

    # Load the reflections from the pickle file
    Command.start('Loading reflections from %s' % args[0])
    table = flex.reflection_table.from_pickle(args[0])
    Command.end('Loaded %d reflections from %s' % (len(table), args[0]))

    # Save the reflections to the HDF5 file
    Command.start('Saving %d reflections to %s' % (len(table), args[1]))
    table.as_h5(args[1])
    Command.end('Saved %d reflections to %s' % (len(table), args[1]))
