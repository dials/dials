#!/usr/bin/env python
#
# dials.process.py
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
        usage = "usage: %prog [options] [param.phil] sweep.json crystal.json"

        # Initialise the base class
        ScriptRunner.__init__(self, usage=usage)

        # Output filename option
        self.config().add_option(
            '--sweep-filename',
            dest = 'sweep_filename',
            type = 'string', default = 'refined_sweep.json',
            help = 'Set the filename for refined experimental models.')

        # Output filename option
        self.config().add_option(
            '--crystal-filename',
            dest = 'crystal_filename',
            type = 'string', default = 'refined_crystal.json',
            help = 'Set the filename for refined crystal model.')

        # Output filename option
        self.config().add_option(
            '--observed-filename',
            dest = 'observed_filename',
            type = 'string', default = 'observed.pickle',
            help = 'Set the filename for observed reflections.')

        # Output filename option
        self.config().add_option(
            '--reference-filename',
            dest = 'reference_filename',
            type = 'string', default = 'reference.pickle',
            help = 'Set the filename for reference profiles.')

        # Output filename option
        self.config().add_option(
            '--integrated-filename',
            dest = 'integrated_filename',
            type = 'string', default = 'integrated.pickle',
            help = 'Set the filename for integrated reflections.')

        # Add a verbosity option
        self.config().add_option(
            "-v", "--verbosity",
            action="count", default=1,
            help="set verbosity level; -vv gives verbosity level 2")


    def main(self, params, options, args):
        '''Execute the script.'''
        from dials.model.serialize import load, dump
        from dials.process import ProcessFactory
        from dials.util.command_line import Command
        from dials.algorithms import shoebox

        # Get the filenames
        if len(args) == 2:
            sweep = load.sweep(args[0])
            crystal = load.crystal(args[1])
        else:
            self.config().print_help()
            return

        # Creating the processing function
        pipeline = ProcessFactory.from_parameters(params, options.verbosity)

        # Process the data
        pipeline.run(sweep, crystal)
        observed = pipeline.observed
        reference = pipeline.reference
        integrated = pipeline.integrated

        # Print some help
        print '\n' + '=' * 80
        print 'Writing results to file'
        print '=' * 80 + '\n'

        # Save the refined geometry to file
        sweep_filename = options.sweep_filename
        Command.start('Saving refined geometry to {0}'.format(sweep_filename))
        dump.sweep(sweep, open(sweep_filename, 'w'))
        Command.end('Saved refined geometry to {0}'.format(sweep_filename))

        # Save the refined crystal to file
        crystal_filename = options.crystal_filename
        Command.start('Saving refined geometry to {0}'.format(crystal_filename))
        dump.crystal(crystal, open(crystal_filename, 'w'))
        Command.end('Saved refined geometry to {0}'.format(crystal_filename))

        # Save the observed reflections to file
        Command.start('Saving observed to {0}'.format(options.observed_filename))
        dump.reflections(observed, options.observed_filename)
        Command.end('Saved observed to {0}'.format(options.observed_filename))

        # Save the reference profiles to file
        Command.start('Saving reference to {0}'.format(options.reference_filename))
        dump.reference(reference, options.reference_filename)
        Command.end('Saved reference to {0}'.format(options.reference_filename))

        # Save the integrated reflections to file
        Command.start('Saving integrated to {0}'.format(options.integrated_filename))
        shoebox.deallocate(integrated)
        dump.reflections(integrated, options.integrated_filename)
        Command.end('Saved integrated to {0}'.format(options.integrated_filename))


if __name__ == '__main__':
    script = Script()
    script.run()
