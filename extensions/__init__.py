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

from __future__ import absolute_import, division


def import_sub_modules(paths):
  '''
  Import all sub modules.

  :param paths: The list of module paths
  :returns: A list of python modules

  '''
  from pkgutil import walk_packages

  # Create a generator to walk through the sub packages
  walker = walk_packages(path=paths, onerror=lambda x: None)

  # Walk through, importing the packages
  modules = []
  for importer, modname, ispkg in walker:
    loader = importer.find_module(modname)
    modules.append(loader.load_module(modname))
  return modules

def import_extensions():
  '''
  Import the extensions

  :returns: The modules containing the python extensions

  '''
  from dials.framework import env
  from os.path import dirname

  # Get the paths to import from
  paths = [dirname(__file__)] + env.cache.paths()

  # Import the sub modules
  return import_sub_modules(paths)

def grab_extensions(modules):
  '''
  Get a list of extensions.

  :param modules: The list of modules
  :returns: The list of extensions

  '''
  from dials.framework.interface import Interface
  from inspect import getmembers, isclass
  ext = []
  for mod in modules:
    for name, obj in getmembers(mod):
      if isclass(obj) and issubclass(obj, Interface):
        ext.append((name, obj))
  return ext


# Import sub modules
_modules = import_extensions()

# Add this extensions to the current namespace
for name, obj in grab_extensions(_modules):
  globals()[name] = obj
