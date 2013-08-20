#!/usr/bin/env python
#
# dials.find_spots.py
#
#  Copyright (C) 2013 Diamond Light Source
#
#  Author: James Parkhurst
#
#  This code is distributed under the BSD license, a copy of which is
#  included in the root directory of this package.
from __future__ import division
from dials.util.script import ScriptRunner


class Script(ScriptRunner):
    '''A class for running the script.'''

    def __init__(self):
        '''Initialise the script.'''

        # The script usage
        usage = "usage: %prog [options] [param.phil] sweep.json"

        # Initialise the base class
        ScriptRunner.__init__(self, usage=usage)

        # Output filename option
        self.config().add_option(
            '-o', '--output-filename',
            dest = 'output_filename',
            type = 'string', default = 'spots.pickle',
            help = 'Set the filename for found spots.')

    def main(self, params, options, args):
        '''Execute the script.'''
        from dials.algorithms.peak_finding.spotfinder_factory \
            import SpotFinderFactory
        from dials.model.serialize import load, dump
        from dials.util.command_line import Command

        # Check the number of arguments is correct
        if len(args) == 1:
            print 'Loading initial models from {0}'.format(args[0])
            sweep = load.sweep(args[0])
        else:
            self.config().print_help()
            return

        # Get the integrator from the input parameters
        print 'Configuring spot finder from input parameters'
        find_spots = SpotFinderFactory.from_parameters(params)

        # Find the strong spots in the sweep
        print 'Finding strong spots'
        reflections = find_spots(sweep)

        # Save the reflections to file
        Command.start('Saving {0} spots to {1}'.format(
            len(reflections), options.output_filename))
        dump.reflections(reflections, options.output_filename)
        Command.end('Saved {0} spots to {1}'.format(
            len(reflections), options.output_filename))


if __name__ == '__main__':
    script = Script()
    script.run()
