#!/usr/bin/env python
#
# options.py
#
#  Copyright (C) 2013 Diamond Light Source
#
#  Author: James Parkhurst
#
#  This code is distributed under the BSD license, a copy of which is
#  included in the root directory of this package.

from __future__ import division
from dials.util import HalError
import optparse


class PhilToDict(object):
  '''Get a dictionary from the phil objects.'''

  def __init__(self):
    '''Init the class'''
    pass

  def __call__(self, phil):
    '''Recursively get the dictionary objects.'''

    d = {}
    for obj in phil.objects:
      if obj.is_scope:
        d[obj.name] = self.__call__(obj)
      else:
        if obj.multiple:
          try:
            d[obj.name] = d[obj.name] + [obj.extract()]
          except Exception:
            d[obj.name] = [obj.extract()]
        else:
          d[obj.name] = obj.extract()

    return d


class OptionParser(optparse.OptionParser):
  '''A class to parse command line options and get the system configuration.
  The class extends optparse.OptionParser to include the reading of phil
  parameters.'''

  def __init__(self, **kwargs):
    '''Initialise the class.'''
    from libtbx.phil import parse

    # Set the system phil scope
    self._system_phil = kwargs.get('phil', parse(""))

    # Set the working phil scope
    self._phil = self._system_phil.fetch(source=parse(""))

    # Delete the phil keyword from the keyword arguments
    del kwargs['phil']

    # Initialise the option parser
    optparse.OptionParser.__init__(self, **kwargs)

    # Add an option to show configuration parameters
    if self._system_phil.as_str() != '':
      self.parser.add_option(
        '-c',
        action='count',
        default=0,
        dest='show_config',
        help='Show the configuration parameters.')

    # Set a verbosity parameter
    self.parser.add_option(
      '-v',
      action='count',
      default=0,
      dest='verbose',
      help='Set the verbosity')

  def parse_args(self, show_diff_phil=False):
    '''Parse the command line arguments and get system configuration.'''
    from dials.framework.config import CommandLineConfig
    # Parse the command line arguments, this will separate out
    # options (e.g. -o, --option) and positional arguments, in
    # which phil options will be included.
    options, args = optparse.OptionParser.parse_args(self)

    # Parse the command line phil parameters
    config = CommandLineConfig(self.system_phil().command_line_argument_interpreter())
    args, command_line_phil = config.parse(args)
    self._phil = self.system_phil().fetch(sources=command_line_phil)
    params = self._phil.extract()

    # Show config
    if options.show_config > 0:
      self.print_phil(attributes_level=options.show_config-1)
      exit(0)

    # Print the diff phil
    if show_diff_phil:
      diff_phil_str = self.diff_phil().as_str()
      if (diff_phil_str is not ''):
        print 'The following parameters have been modified:\n'
        print diff_phil_str

    # Return the parameters
    return params, options, args

  def phil(self):
    '''Get the phil object'''
    return self._phil

  def system_phil(self):
    '''Get the system phil.'''
    return self._system_phil

  def diff_phil(self):
    ''' Get the diff phil. '''
    return self.system_phil().fetch_diff(source=self.phil())

  def print_phil(self, attributes_level=0):
    '''Print the system and command line parameters.'''
    print self.phil().as_str(attributes_level=attributes_level)

  def print_system_phil(self, attributes_level=0):
    '''Print the system parameters.'''
    print self.system_phil().as_str(attributes_level=attributes_level)

  def print_diff_phil(self, attributes_level=0):
    ''' Print the diff parameters. '''
    print self.diff_phil().as_str(attributes_level=attributes_level)

  def format_epilog(self, formatter):
    ''' Don't do formatting on epilog. '''
    if self.epilog is None:
      return ''
    return self.epilog

