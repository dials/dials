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

class ScriptRunner(object):
  '''Run the script.'''

  def __init__(self, **kwargs):
    '''Set the input and output filenames.'''

    self.pickle_filename = kwargs['pickle_filename']
    self.nexus_filename = kwargs['nexus_filename']

  def run(self):
    '''Load reflections from pickle file and save to HDF5 file.'''

    from dials.model.data import Reflection, ReflectionList
    from dials.util.nexus import NexusFile
    import cPickle as pickle

    # Load the reflections from the pickle file
    print 'Loading reflections from {0}'.format(self.pickle_filename)
    reflections = pickle.load(open(self.pickle_filename, 'rb'))

    print 'Saving {0} reflection to {1}'.format(
        len(reflections), self.nexus_filename)
    nexus = NexusFile(self.nexus_filename, 'w')
    nexus.set_reflections(reflections)
    nexus.close()


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
