
class CompletionGenerator(object):
    '''A class to auto-generate bash-completion scripts for dials'''

    def __init__(self):
        '''Initialise the class with the list of programs'''
        import libtbx.load_env
        import os

        # A dictionary of program/phil scope pairs
        self.programs = {
            'spotfinder' : 'spotfinder',
            'subtract_background' : 'background'
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
            text += '\n\n' + self._generate_single(p, options[s])

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
            options[scope].append(option)

        # Return the options
        return options

    def _generate_single(self, program, options):
        '''Generate a single script from the template.'''

        from string import Template

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
