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

    # Try to get the home scope
    self._system_phil = kwargs['phil']
    self._phil = self._system_phil.fetch(source=parse(""))
    del(kwargs['phil'])

    # Initialise the option parser
    optparse.OptionParser.__init__(self, **kwargs)

  def parse_args(self):
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
    return self.epilog
