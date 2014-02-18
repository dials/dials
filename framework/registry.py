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

  def __init__(self):

    # Load the interfaces and extensions
    from dials.framework.interface import Interface
    import dials.interfaces
    import dials.extensions
    self._interface = Interface

    # Read the system config
    from dials.util.options import SystemConfig
    self._config = SystemConfig()
    self._params = self._config.config().extract()

  def __iter__(self):
    return self.interfaces()

  def __len__(self):
    return len(list(self.interfaces()))

  def __contains__(self, iface):
    return iface in list(self)

  def params(self):
    return self._params

  def interfaces(self):
    return self._interface.interfaces()

  def extensions(self, iface):
    return iface.extensions()

  def all_extensions(self):
    for iface in self:
      for ext in self.extensions(iface):
        yield ext

  def interface(self, name):
    return dict((iface.name, iface) for iface in self.interfaces())[name]

  def __getitem__(self, key):
    name = getattr(self._params, key).algorithm
    return self.interface(key).extension(name)


def init_ext(iface, *args, **kwargs):
  ''' Helper function to instantiate an extension. '''
  registry = Registry()
  return registry[iface](registry.params(), *args, **kwargs)

