#!/usr/bin/env python
#
# dials.reference_profiles.py
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
        usage += "sweep.json crystal.json strong.pickle"

        # Initialise the base class
        ScriptRunner.__init__(self, usage=usage)

        # Output filename option
        self.config().add_option(
            '-o', '--output-filename',
            dest = 'output_filename',
            type = 'string', default = 'reference.pickle',
            help = 'Set the filename for the reference profiles.')

    def main(self, params, options, args):
        '''Execute the script.'''
        from dials.algorithms.integration import ReferenceProfileFactory
        from dials.model.serialize import load, dump
        from dials.util.command_line import Command

        # Check the number of arguments is correct
        if len(args) != 3:
            self.config().print_help()
            return

        # Load the sweep
        Command.start('Loading sweep from {0}'.format(args[0]))
        sweep = load.sweep(args[0])
        Command.end('Loaded sweep from {0}'.format(args[0]))

        # Load the crystal
        Command.start('Loading crystal from {0}'.format(args[1]))
        crystal = load.crystal(args[1])
        Command.end('Loaded crystal from {0}'.format(args[1]))

        # Load the strong reflections
        Command.start('Loading strong spots from {0}'.format(args[2]))
        strong = load.reflections(args[2])
        Command.end('Loaded {0} strong spots from {1}'.format(
            len(strong), args[2]))

        # Get the reference profile creator from the input parameters
        print 'Configuring reference profile creator from input parameters'
        create_reference = ReferenceProfileFactory.from_parameters(params)

        # Create the reference profiles
        print 'Create the reference profiles'
        reference = create_reference(sweep, crystal, strong)

        # Save the reference profiles to file
        Command.start('Saving {0} reference profiles to {1}'.format(
            len(reference), options.output_filename))
        dump.reference(reference, options.output_filename)
        Command.end('Saved {0} reference profiles to {1}'.format(
            len(reference), options.output_filename))


if __name__ == '__main__':
    script = Script()
    script.run()
