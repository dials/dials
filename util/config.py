#!/usr/bin/env python
#
# config.py
#
#  Copyright (C) 2013 Diamond Light Source
#
#  Author: James Parkhurst
#
#  This code is distributed under the BSD license, a copy of which is
#  included in the root directory of this package.

from __future__ import division

class CompletionGenerator(object):
    '''A class to auto-generate bash-completion scripts for dials'''

    def __init__(self):
        '''Initialise the class with the list of programs'''
        import libtbx.load_env
        import os

        # A dictionary of program/phil scope pairs
        self.programs = {
            'parameters' : None,
            'spotfinder' : ['spotfinder', 'logging', 'lookup'],
            'integrate' : ['integration', 'logging', 'lookup'],
            'refine' : ['refinement', 'logging'],
        }

        # Find the dials distribution directory
        dist_path = libtbx.env.dist_path('dials')
        build_path = abs(libtbx.env.build_path)

        # Get the location of the template file
        self.template_filename = os.path.join(
            dist_path, 'data', 'completion_template.txt')
        self.output_filename = os.path.join(
            build_path, 'dials', 'completion')

    def generate(self):
        '''Generate the script.'''

        # Get scoped options
        options = self._get_scoped_options()

        # Generate a script for each program
        text = '#!/bin/bash'
        for p, s in self.programs.iteritems():
            if s:
                opt = [o for ol in s for o in options[ol]]
            else:
                opt = [o for ol in options for o in options[ol]]
            text += '\n\n' + self._generate_single(p, opt)

        # Save the generated completion stuff to file
        with open(self.output_filename, 'w') as f:
            f.write(text)

    def _get_scoped_options(self):
        '''Scope the phil options'''

        from collections import defaultdict
        from dials.util.options import SystemConfig

        # Get the system configuration
        system = SystemConfig()
        options = system.config()

        # Get all the definitions
        definitions = options.all_definitions()

        # Get all the options per top-level scope
        options = defaultdict(list)
        for d in definitions:
            path = d.path
            scope, option = path.split('.', 1)
            options[scope].append(path)

        # Return the options
        return options

    def _generate_single(self, program, options):
        '''Generate a single script from the template.'''

        from string import Template

        # Check the number of options
        if len(options) == 0:
            return ''

        # Open the completion template
        with open(self.template_filename, 'r') as f:
            text = f.read()

        # Create the template
        template = Template(text)

        # Indent in template for options
        indent = 7

        # Create the option string
        option_str = reduce(lambda a, b: a + '\n' + ' ' * indent + b, options)

        # Substitute the variables
        return template.substitute(program=program, options=option_str)


if __name__ == '__main__':
    gen = CompletionGenerator()
    gen.generate()
