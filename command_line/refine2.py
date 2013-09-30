#!/usr/bin/env python
#
# dials.refine.py
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
        usage  = "usage: %prog [options] [param.phil] " \
                 "sweep.json crystal.json reflections.pickle"

        # Initialise the base class
        ScriptRunner.__init__(self, usage=usage)

        # Output filename option
        self.config().add_option(
            '--output-sweep-filename',
            dest = 'output_sweep_filename',
            type = 'string', default = 'refined_sweep.json',
            help = 'Set the filename for refined experimental models.')

        # Output filename option
        self.config().add_option(
            '--output-crystal-filename',
            dest = 'output_crystal_filename',
            type = 'string', default = 'refined_crystal.json',
            help = 'Set the filename for refined crystal model.')

        # Add a verbosity option
        self.config().add_option(
            "-v", "--verbosity",
            action="count", default=0,
            help="set verbosity level; -vv gives verbosity level 2")

    def main(self, params, options, args):
        '''Execute the script.'''
        from dials.algorithms.refinement import RefinerFactory
        from dials.model.serialize import load, dump
        from dials.model.data import ReflectionList
        import cPickle as pickle

        # Check the number of arguments is correct
        if len(args) != 3:
            self.config().print_help()
            return

        # Get the refiner from the input parameters
        print 'Configurating refiner from input parameters'
        refine = RefinerFactory.from_parameters(params, options.verbosity)

        # Try to load the models
        print 'Loading models from {0} and {1}'.format(args[0], args[1])
        sweep = load.sweep(args[0])
        crystal = load.crystal(args[1])
        reflections = pickle.load(open(args[2], 'rb'))

        # Refine the geometry
        print 'Performing refinement'
        refined = refine(sweep, crystal, reflections)

        # The sweep and crystal are updated by side-effect

        # Save the refined geometry to file
        output_sweep_filename = options.output_sweep_filename
        print 'Saving refined geometry to {0}'.format(output_sweep_filename)
        dump.sweep(sweep, open(output_sweep_filename, 'w'))

        # Save the refined crystal to file
        output_crystal_filename = options.output_crystal_filename
        print 'Saving refined geometry to {0}'.format(output_crystal_filename)
        dump.crystal(crystal, open(output_crystal_filename, 'w'))


if __name__ == '__main__':
    script = Script()
    script.run()
