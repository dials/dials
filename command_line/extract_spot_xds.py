#!/usr/bin/env python
#
# extract_spot_xds.py
#
#  Copyright (C) 2013 Diamond Light Source
#
#  Author: James Parkhurst
#
#  This code is distributed under the BSD license, a copy of which is
#  included in the root directory of this package.


class ScriptRunner(object):
    '''Class to run script.'''

    def __init__(self, spot_filename, output_filename, include_invalid):
        '''Setup the script.'''

        # Filename data
        self.spot_filename = spot_filename
        self.output_filename = output_filename
        self.include_invalid = include_invalid

    def __call__(self):
        '''Run the script.'''
        from iotbx.xds import spot_xds
        from dials.model.data import Reflection, ReflectionList
        from dials.util.command_line import Command

        # Read the SPOT.XDS file
        Command.start('Reading SPOT.XDS')
        handle = spot_xds.reader()
        handle.read_file(self.spot_filename)

        # Get the info
        centroid = handle.centroid
        intensity = handle.intensity
        try:
            miller_index = handle.miller_index
        except AttributeError:
            miller_index = None

        Command.end('Read {0} spots from SPOT.XDS file.'.format(len(centroid)))

        # Create the reflection list
        Command.start('Creating reflection list')
        if miller_index:
            rlist = ReflectionList(len(centroid))
            for r, c, i, h in zip(rlist, centroid, intensity, miller_index):
                if not self.include_invalid and h != (0, 0, 0):
                    r.centroid_position = c
                    r.centroid_variance = (1.0, 1.0, 1.0)
                    r.intensity = i
                    r.miller_index = h

        else:
            rlist = ReflectionList(len(centroid))
            for r, c, i in zip(rlist, centroid, intensity):
                r.centroid_position = c
                r.centroid_variance = (1.0, 1.0, 1.0)
                r.intensity = i
        Command.end('Created reflection list')

        # Save the reflection list
        if self.output_filename != None:
            import pickle
            Command.start('Saving reflections to {0}'.format(
                self.output_filename))
            pickle.dump(rlist, open(self.output_filename, 'wb'))
            Command.end('Saved reflections to {0}'.format(
                self.output_filename))


if __name__ == '__main__':

    from optparse import OptionParser

    # Specify the command line options
    usage  = "usage: %prog [options] /path/to/SPOT.XDS"

    # Create an option parser
    parser = OptionParser(usage)

    # Add command line options
    parser.add_option('-z', '--zero-hkl',
                      dest='zero_hkl', action="store_true", default=False,
                      help='Include invalid reflections')
    parser.add_option('-o', '--output-file',
                      dest='output_file', type="string", default="",
                      help='Destination filename for reflections')

    # Parse the arguments
    options, args = parser.parse_args()

    # Print help if no arguments specified, otherwise call function
    if len(args) < 1:
        parser.print_help()
    else:
        # Initialise the script runner
        runner = ScriptRunner(
            spot_filename=args[0],
            output_filename=options.output_file,
            include_invalid=options.zero_hkl)

        # Run the script
        runner()
