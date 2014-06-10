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


class CommandLineConfig(object):
  '''A class to read the commandline configuration.'''

  def __init__(self, interpretor):
    '''Set the command line interpretor.'''
    self._interpretor = interpretor

  def parse(self, argv):
    '''Get the configuration.'''
    import os

    phils, positionals = self._interpretor.process_and_fetch(
      argv, custom_processor="collect_remaining")
    # Return positional arguments and phils
    return positionals, [phils]


class ConfigWriter(object):
  '''Class to write configuration to file.'''

  def __init__(self, master_phil):
    '''Initialise with the master phil.'''
    self._master_phil = master_phil

  def write(self, params, filename):
    '''Write the configuration to file.'''
    # Get the modified phil
    modified_phil = self._master_phil.format(python_object=params)

    # Get the phil as text
    text = modified_phil.as_str()

    # Write the text to file
    with open(filename, 'w') as f:
      f.write(text)


class Config(object):
  ''' Manage configuration. '''

  def __init__(self):
    ''' Initialise by reading the system phil. '''
    from libtbx.phil import parse

    # Parse all the phil files
    self._system_phil = parse(
    '''
      include scope dials.data.spotfinding.phil_scope
      include scope dials.data.integration.phil_scope
      include scope dials.data.refinement.phil_scope
      include scope dials.data.indexing.phil_scope
    ''', process_includes=True)

    self._phil = self._system_phil
    self._params = self._phil.extract()

  def system_phil(self):
    ''' Return the system phil. '''
    return self._system_phil

  def phil(self):
    ''' Return the user phil. '''
    return self._phil

  def params(self):
    ''' Return the cached parameters. '''
    return self._params

  def write(self, filename):
    ''' Write the configuration to file. '''
    writer = ConfigWriter(self._system_phil)
    writer.write(self.params(), filename)

  def fetch(self, sources=None, source=None):
    ''' Do the phil fetch. '''
    self._phil = self._system_phil.fetch(source=source, sources=sources)
    self._params = self._phil.extract()

  def parse(self, text):
    ''' Parse the given phil string. '''
    from iotbx.phil import parse
    self._phil = self.fetch(source=parse(text))

  def try_parse(self, args):
    ''' Try to parse arguments and return unhandled. '''
    config = CommandLineConfig(
      self._system_phil.command_line_argument_interpreter())
    args, phil = config.parse(args)
    self.fetch(sources=phil)
    return args
