#!/usr/bin/env python
#
# display_spots.py
#
#  Copyright (C) 2013 Diamond Light Source
#
#  Author: Richard Gildea
#
#  This code is distributed under the BSD license, a copy of which is
#  included in the root directory of this package.
from __future__ import division
# LIBTBX_PRE_DISPATCHER_INCLUDE_SH export PHENIX_GUI_ENVIRONMENT=1
# LIBTBX_PRE_DISPATCHER_INCLUDE_SH export BOOST_ADAPTBX_FPE_DEFAULT=1


class ScriptRunner(object):
    '''Class to run script.'''

    def __init__(self, sweep_filenames, reflections_filename):
        '''Setup the script.'''

        # Filename data
        self.sweep_filenames = sweep_filenames
        self.reflections_filename = reflections_filename

    def __call__(self):
        '''Run the script.'''
        import cPickle as pickle
        from dials.model.data import ReflectionList # import dependency
        from dials.util.command_line import Command

        # Read the pickle file
        Command.start('Reading reflection file.')
        with open(self.reflections_filename, 'rb') as f:
            self.reflections = pickle.load(f)

        Command.end('Read {0} spots from reflection file.'.format(
            len(self.reflections)))

        self.view()

    def view(self):
        from dials.util.spotfinder_wrap import spot_wrapper
        spot_wrapper(working_phil=None).display(
            sweep_filenames=self.sweep_filenames, reflections=self.reflections)

if __name__ == '__main__':

    from optparse import OptionParser

    # Specify the command line options
    usage  = "usage: %prog [options] "
    usage += "/path/to/reflections.pickle "
    usage += "/path/to/image/files "

    # Create an option parser
    parser = OptionParser(usage)

    # Parse the arguments
    options, args = parser.parse_args()

    # Print help if no arguments specified, otherwise call function
    if len(args) < 1:
        parser.print_help()
    else:
        # Initialise the script runner
        runner = ScriptRunner(
            reflections_filename=args[0],
            sweep_filenames=args[1:])

        # Run the script
        runner()
