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

class LoggingConfig(object):
    '''Class to configure logging.'''

    def __init__(self):
        '''Init the class'''
        pass

    def configure(self, phil):
        '''Configure the logging.'''
        import logging.config
        from dials.util.options import PhilToDict

        # Get the dictionary from the phil spec
        logging_dict = PhilToDict()(phil)['logging']

        # Create the filename and logging name
        self.create_filename(logging_dict)
        self.create_loggername(logging_dict)

        # Configure the logger
        logging.config.dictConfig(logging_dict)

    def create_filename(self, logging_dict):
        '''Hack to set the log file name'''
        from datetime import datetime
        import os

        # Hack to create the file name from the config. Get and create dir
        directory = logging_dict['handlers']['file']['directory']
        directory = os.path.expanduser(directory)
        if not os.path.exists(directory):
            os.makedirs(directory)

        # Get filename
        filename = 'dials_{0}.log'.format(
            datetime.now().strftime("%Y-%m-%d_%H_%M_%S"))
        filename = os.path.expanduser(os.path.join(directory, filename))
        del(logging_dict['handlers']['file']['directory'])
        logging_dict['handlers']['file']['filename'] = filename

    def create_loggername(self, logging_dict):
        '''Hack to change logger name.'''
        logging_dict['loggers'][''] = logging_dict['loggers']['default']
        del(logging_dict['loggers']['default'])


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
        '''Run the script. This function calls the 'main' method that
        should be overloaded in your script.'''
        from dials.util.options import ConfigWriter
        import sys

        # Parse the configuration parameters
        params, options, args = self._config.parse_args()

        # Configure the logging
        self._configure_logging(self._config.phil())

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

    def _configure_logging(self, phil):
        '''Configure the logging.'''
        log = LoggingConfig()
        log.configure(phil)

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
