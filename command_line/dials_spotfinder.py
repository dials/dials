#!/usr/bin/env python
#
# dials.dials_spotfinder.py
#
#  Copyright (C) 2013 Diamond Light Source
#
#  Author: James Parkhurst
#
#  This code is distributed under the BSD license, a copy of which is
#  included in the root directory of this package.

from __future__ import division
from dials.util.script import ScriptRunner
# LIBTBX_SET_DISPATCHER_NAME dials.spotfinder

class Script(ScriptRunner):
    '''A class for running the script.'''

    def __init__(self):
        '''Initialise the script.'''

        # The script usage
        usage = "usage: %prog [options] [param.phil] "\
                "{sweep.json | image1.file [image2.file ...]}"

        # Initialise the base class
        ScriptRunner.__init__(self, usage=usage)

        # Output filename option
        self.config().add_option(
            '-o', '--output-filename',
            dest = 'output_filename',
            type = 'string', default = 'strong.pickle',
            help = 'Set the filename for found strong spots.')

    def main(self, params, options, args):
        '''Execute the script.'''
        from dials.algorithms.peak_finding.spotfinder_factory \
            import SpotFinderFactory
        from dials.algorithms import shoebox
        from dials.model.serialize import load, dump
        from dials.util.command_line import Command
        from dxtbx.imageset import ImageSetFactory
        from dials.util.command_line import Importer

        # Try importing the command line arguments
        importer = Importer(args)
        if len(importer.imagesets) == 0:
            self.config().print_help()
            return
        elif len(importer.imagesets) > 1:
            raise RuntimeError("Only one imageset can be processed at a time")
        sweep = importer.imagesets[0]

        # Get the integrator from the input parameters
        print 'Configuring spot finder from input parameters'
        find_spots = SpotFinderFactory.from_parameters(params)

        # Find the strong spots in the sweep
        print 'Finding strong spots'
        reflections = find_spots(sweep)

        # Save the reflections to file
        Command.start('Saving {0} reflections to {1}'.format(
            len(reflections), options.output_filename))
        dump.reflections(reflections, options.output_filename)
        Command.end('Saved {0} reflections to {1}'.format(
            len(reflections), options.output_filename))


if __name__ == '__main__':
    script = Script()
    script.run()
