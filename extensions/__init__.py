#!/usr/bin/env python
#
# __init__.py
#
#  Copyright (C) 2013 Diamond Light Source
#
#  Author: James Parkhurst
#
#  This code is distributed under the BSD license, a copy of which is
#  included in the root directory of this package.

from __future__ import division


def import_sub_modules():
  ''' Import all sub modules. '''
  from pkgutil import walk_packages
  from os.path import dirname

  # Create a generator to walk through the sub packages
  path = dirname(__file__)
  walker = walk_packages(path=[path], onerror=lambda x: None)

  # Walk through, importing the packages
  modules = []
  for importer, modname, ispkg in walker:
    modules.append(__import__(modname, globals()))
  return modules

def grab_extensions(modules):
  ''' Get a list of extensions. '''
  from dials.framework.interface import Interface
  from inspect import getmembers, isclass
  ext = []
  for mod in modules:
    for name, obj in getmembers(mod):
      if isclass(obj) and issubclass(obj, Interface):
        ext.append((name, obj))
  return ext


# Import sub modules
_modules = import_sub_modules()

# Add this extensions to the current namespace
for name, obj in grab_extensions(_modules):
  globals()[name] = obj
