#!/usr/bin/env python
#
# plugin.py
#
#  Copyright (C) 2013 Diamond Light Source
#
#  Author: James Parkhurst
#
#  This code is distributed under the BSD license, a copy of which is
#  included in the root directory of this package.
from __future__ import division


def singleton(cls):
  ''' Make a class a singleton. '''
  import functools
  instance = cls()
  @functools.wraps(cls)
  def wrapper():
    return instance
  return wrapper


@singleton
class Registry(object):
  ''' Registry to help loading extensions. '''

  def __init__(self):

    # Load the interfaces and extensions
    from dials.framework.interface import Interface
    import dials.interfaces
    import dials.extensions
    self._interface = Interface

    # Read the system config
    from dials.framework.config import Config
    self._config = Config()

  def __iter__(self):
    ''' Iterate through the interfaces. '''
    return self.interfaces()

  def __len__(self):
    ''' Get the number of interfaces. '''
    return len(list(self.interfaces()))

  def __contains__(self, iface):
    ''' Check the registry contains the interface. '''
    return iface in list(self)

  def config(self):
    ''' Get the configuration. '''
    return self._config

  def params(self):
    ''' Get the parameters. '''
    return self._config.params()

  def interfaces(self):
    ''' Get the interfaces. '''
    return self._interface.interfaces()

  def extensions(self, iface):
    ''' Get the extensions for an interface. '''
    return iface.extensions()

  def all_extensions(self):
    ''' Get a list of all extensions. '''
    for iface in self:
      for ext in self.extensions(iface):
        yield ext

  def interface(self, name):
    ''' Get an interface. '''
    return dict((iface.name, iface) for iface in self.interfaces())[name]

  def __getitem__(self, key):
    ''' Get an algorithm. '''
    path = key.split('.')
    obj = self.params()
    for p in path:
      obj = getattr(obj, p)
    name = obj.algorithm
    return self.interface(key).extension(name)


def init_ext(iface, *args, **kwargs):
  ''' Helper function to instantiate an extension. '''
  registry = Registry()
  return registry[iface](registry.params(), *args, **kwargs)

