#!/usr/bin/env python
#
# dials.util.script.py
#
#  Copyright (C) 2013 Diamond Light Source
#
#  Author: James Parkhurst
#
#  This code is distributed under the BSD license, a copy of which is
#  included in the root directory of this package.

class ScriptRunner(object):
    '''A class to write scripts in a consistent way.'''

    def __init__(self, **kwargs):
        '''Initialise the class.'''

        # Create the options
        self._create_options(**kwargs)

    def _create_options(self, **kwargs):
        '''Create the generic options for the script.'''
        from dials.util.options import OptionParser

        # Create the option parser
        self._config = OptionParser(
            usage      = kwargs.get('usage', ''),
            home_scope = kwargs.get('home_scope', ''))

        # Add an option to show the current configuraton
        if kwargs.get('show_config_option', True):
            self._config.add_option(
                '--show-config',
                dest='show_config', action='store_true', default=False,
                help='Show the configuration parameters.')

        # Add an option to save the working configuration
        if kwargs.get('save_config_option', True):
            self._config.add_option(
                '--save-config',
                dest='output_config_file', type='string', default=None,
                help='The file in which to save the configuration parameters.')

        # Bind the config add_option function to this class
        self.add_option = self._config.add_option

    def run(self):
        '''Run the script. This function calls the 'start' method that
        should be overloaded in your script.'''
        from dials.util.options import ConfigWriter
        import sys

        # Parse the configuration parameters
        params, options, args = self._config.parse_args()

        # Save the working parameters
        self.working_params = params

        # If the show_config flag is set, then print and exit
        try:
            if options.show_config:
                self._config.print_phil()
                sys.exit(0)
        except AttributeError:
            pass

        # Run the actual script
        self.main(params, options, args)

        # If the save_config flag is set, then save the configuration
        try:
            if options.output_config_file:
                writer = ConfigWriter(self._config.system_phil())
                writer.write(self.working_params, options.output_config_file)
        except AttributeError:
            pass

    def config(self):
        '''Get the configuration'''
        return self._config

    def update_params(self, params):
        '''Update the working parameters'''
        self.working_params=params

    def main(self, params, options, args):
        '''The main body of the script: overload!

        Params:
            params The phil parameters
            options The command line options
            args The command line positional arguments

        '''
        raise RuntimeError('Overload!')


if __name__ == '__main__':

    script = ScriptRunner()
    script.run()
