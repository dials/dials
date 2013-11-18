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

class SystemConfigReader(object):
  '''Class to read system configuration.'''

  def __init__(self):
    '''Initialise the class.'''
    pass

  def master(self):
    '''Get the master config text.'''
    return '\n'.join([self._read_file(filename, True)
        for filename in self.master_filenames()])

  def user(self):
    '''Get the user config text.'''
    return self._read_file(self.user_filename(), False)

  def _read_file(self, filename, fail=False):
    '''Read the config file.'''

    # Try to read the phil file
    text = ''
    if filename:
      try:
        with open(filename, 'r') as f:
          text = f.read()
      except IOError:
        if fail:
          raise HalError('error reading {0}'.format(filename))

    # Return the text
    return text

  def master_filenames(self):
    '''Get the master filename.'''
    import libtbx.load_env
    import os

    # Find the dials distribution directory
    try:
      path = libtbx.env.dist_path('dials')
    except KeyError, e:
      raise HalError('dials is not configured.')

    # Get the location of the master file
    return [os.path.join(path, 'data', 'logging.phil'),
            os.path.join(path, 'data', 'lookup.phil'),
            os.path.join(path, 'data', 'spotfinding.phil'),
            os.path.join(path, 'data', 'integration.phil'),
            os.path.join(path, 'data', 'refinement.phil')]

  def user_filename(self):
    '''Get the user filename.'''
    import os

    # Get the location of the user filename
    return os.path.join(os.path.expanduser('~'), '.dialsrc')


class SystemConfig(object):
  '''A class to read the system configuration.'''

  def __init__(self):
    '''Initialise the class.'''

    # Create the config file reader
    self._files = SystemConfigReader()

  def config(self):
    '''Get the configuration.'''
    from libtbx.phil import parse

    # Read the master and user files
    master_text = self._files.master()
    user_text   = self._files.user()

    # Parse the files with phil
    try:
      self.master_phil  = parse(master_text)
    except RuntimeError, e:
      raise HalError(str(e) + ' ({0})'.format(
          self._files.master_filenames()))

    # Parse the user file with phil
    try:
      self.user_phil = parse(user_text)
    except RuntimeError, e:
      raise HalError(str(e) + ' ({0})'.format(
          self._files.user_filename()))

    # Fetch the working phil from all the system sources
    try:
      phil = self.master_phil.fetch(sources = [self.user_phil])
    except Exception, e:
      raise HalError(e)

    # Return the system phil
    return phil


class CommandLineConfig(object):
  '''A class to read the commandline configuration.'''

  def __init__(self, interpretor):
    '''Set the command line interpretor.'''
    self._interpretor = interpretor

  def config(self, argv):
    '''Get the configuration.'''
    import os

    # For each command line option, try to parse using phil, if an
    # exception occurs, then append to positional arguments, otherwise
    # append to phils
    positionals, phils = [], []
    for arg in argv:
      if os.path.isfile(arg):
        try:
          phils.append(file_name=self._interpretor.process_arg(arg))
        except Exception:
          positionals.append(arg)
      elif arg.find('=') > 0:
        try:
          phils.append(self._interpretor.process_arg(arg))
        except Exception, e:
          raise HalError(e)

    # Return positional arguments and phils
    return positionals, phils


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

    # Try to get the home scope
    try:
      self._scope = kwargs['home_scope']
      del(kwargs['home_scope'])
    except KeyError:
      self._scope = ''

    # Initialise the option parser
    optparse.OptionParser.__init__(self, **kwargs)

  def parse_args(self):
    '''Parse the command line arguments.'''
    import sys

    # Parse arguments, if an error occurs then exit
    try:
      result = self._parse_and_extract_args()
    except HalError, e:
      print e
      sys.exit(0)

    # Return the result
    return result

  def _parse_and_extract_args(self):
    '''Parse and extract the arguments.'''

    # Parse the arguments
    phil, options, args = self._parse_args_internal()

    # Extract the parameters
    try:
      parameters = phil.extract()
    except RuntimeError, e:
      raise HalError(e)

    # return the result
    return parameters, options, args

  def _parse_args_internal(self):
    '''Parse the command line arguments and get system configuration.'''

    # Parse the command line arguments, this will separate out
    # options (e.g. -o, --option) and positional arguments, in
    # which phil options will be included.
    options, args = optparse.OptionParser.parse_args(self)

    # Parse the system arguments from the master and user files
    sysconfig = SystemConfig()
    self._system_phil = sysconfig.config()

    # Parse the command line arguments. This will separate true
    # positional arguments from phil arguments.
    comconfig = CommandLineConfig(
        self._system_phil.command_line_argument_interpreter(
            home_scope=self._scope))

    args, command_phil = comconfig.config(args)

    # Fetch the working phil from all the sources
    try:
      self._phil = self._system_phil.fetch(sources = command_phil)
    except Exception, e:
      raise HalError(e)

    # Return the parameters
    return self._phil, options, args

  def phil(self):
    '''Get the phil object'''
    return self._phil

  def system_phil(self):
    '''Get the system phil.'''
    return self._system_phil

  def print_phil(self, attributes_level=0):
    '''Print the system and command line parameters.'''
    import sys

    # Parse the system and command line arguments
    try:
      phil, options, args = self._parse_args_internal()
    except HalError, e:
      print e
      sys.exit(0)

    # Print the phil parameters
    print phil.as_str(attributes_level=attributes_level)

  def print_system_phil(self, attributes_level=0):
    '''Print the system parameters.'''
    import sys

    # Parse the system arguments from the master and user files
    try:
      sysconfig = SystemConfig()
      system_phil = sysconfig.config()
    except HalError, e:
      print e
      sys.exit(0)

    # Print the system parameters
    print system_phil.as_str(attributes_level=attributes_level)


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


if __name__ == '__main__':

  parser = OptionParser(home_scope='spotfinder')
  parser.add_option('-a', dest='a', type='int', help='hello world')
  print parser.parse_args()
  parser.print_phil(attributes_level=1)
  parser.print_system_phil()
