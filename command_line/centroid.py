#!/usr/bin/env python
#
# dials.centroid.py
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
        usage  = "usage: %prog [options] [param.phil] "
        usage += "sweep.json crystal.json [predicted.pickle]"

        # Initialise the base class
        ScriptRunner.__init__(self, usage=usage)

        # Output filename option
        self.config().add_option(
            '-o', '--output-filename',
            dest = 'output_filename',
            type = 'string', default = 'centroided.pickle',
            help = 'Set the filename for centroided reflections.')


    def main(self, params, options, args):
        '''Execute the script.'''
        from dials.algorithms.centroid import CentroidFactory
        from dials.algorithms import shoebox
        from dials.model.serialize import load, dump

        # Check the number of arguments is correct
        if len(args) < 2 or len(args) > 3:
            self.config().print_help()
            return

        # Get the integrator from the input parameters
        print 'Configurating centroid from input parameters'
        centroid = CentroidFactory.from_parameters(params)

        # Try to load the models
        print 'Loading models from {0} and {1}'.format(args[0], args[1])
        sweep = load.sweep(args[0])
        crystal = load.crystal(args[1])

        # Try to load the predicted reflections
        if len(args) == 3:
            print 'Loading predictions from {0}'.format(args[2])
            predicted = load.reflections(args[2])
        else:
            predicted = None

        # Intregate the sweep's reflections
        print 'Finding strong spot centroids'
        reflections = centroid(sweep, crystal, predicted)

        # Save the reflections to file
        print 'Saving reflections to {0}'.format(options.output_filename)
        dump.reflections(reflections, options.output_filename)


if __name__ == '__main__':
    script = Script()
    script.run()
