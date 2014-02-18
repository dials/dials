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
  for importer, modname, ispkg in walker:
    __import__(modname, globals())

# Import sub modules
import_sub_modules()
