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

    # Try to get the home scope
    try:
      self._scope = kwargs['home_scope']
      del(kwargs['home_scope'])
    except KeyError:
      self._scope = ''

    # Initialise the option parser
    optparse.OptionParser.__init__(self, **kwargs)

  def parse_args(self):
    '''Parse the command line arguments and get system configuration.'''
    from dials.framework.registry import Registry
    registry = Registry()

    # Parse the command line arguments, this will separate out
    # options (e.g. -o, --option) and positional arguments, in
    # which phil options will be included.
    options, args = optparse.OptionParser.parse_args(self)
    args = registry.config().try_parse(args)
    return registry.params(), options, args

  def phil(self):
    '''Get the phil object'''
    from dials.framework.registry import Registry
    return Registry().config().phil()

  def system_phil(self):
    '''Get the system phil.'''
    from dials.framework.registry import Registry
    return Registry().config().system_phil()

  def print_phil(self, attributes_level=0):
    '''Print the system and command line parameters.'''
    import sys
    from dials.framework.registry import Registry
    print Registry().config().phil().as_str(attributes_level=attributes_level)

  def print_system_phil(self, attributes_level=0):
    '''Print the system parameters.'''
    import sys
    from dials.framework.registry import Registry
    print registry.config().system_phil().as_str(attributes_level=attributes_level)




if __name__ == '__main__':

  parser = OptionParser(home_scope='spotfinder')
  parser.add_option('-a', dest='a', type='int', help='hello world')
  print parser.parse_args()
  parser.print_phil(attributes_level=1)
  parser.print_system_phil()
