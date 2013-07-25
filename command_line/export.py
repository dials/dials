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

    from optparse import OptionParser

    # Specify the command line options
    usage  = "usage: %prog [options] "
    usage += "/path/to/reflection/file.pkl "
    usage += "/path/to/reflection/file.h5"

    # Create an option parser
    parser = OptionParser(usage)

    # Parse the arguments
    options, args = parser.parse_args()

    if len(args) < 2:
        parser.print_help()

    else:

        # Run the script
        script = ScriptRunner(
            pickle_filename=args[0],
            nexus_filename=args[1])
        script.run()
