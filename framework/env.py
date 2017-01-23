#!/usr/bin/env python
#
# env.py
#
#  Copyright (C) 2013 Diamond Light Source
#
#  Author: James Parkhurst
#
#  This code is distributed under the BSD license, a copy of which is
#  included in the root directory of this package.

from __future__ import absolute_import, division


class Cache(object):
  ''' Simple class to maintain a cache of additional folders to check. '''

  def __init__(self):
    ''' Initialise an open the cache for writing. '''
    import libtbx.load_env
    from os.path import join

    # Set the paths
    path = libtbx.env.build_path.sh_value()
    self._path = join(path, "dials_extension_folders")

    # Try to read some paths
    try:
      with open(self._path, "r") as infile:
        self._cache = set(infile.readlines())
    except IOError:
        self._cache = set()

  def add(self, path):
    ''' Add a folder. '''
    print ' Adding %s to dials extension folders' % path
    self._cache.add(path)

  def wipe(self):
    ''' Wipe the cache. '''
    self._cache.clear()

  def paths(self):
    ''' Read the paths. '''
    return list(self._cache)

  def __del__(self):
    ''' Write on out. '''
    try:
      with open(self._path, "w") as outfile:
        outfile.write('\n'.join(self._cache))
    except Exception:
      pass


# Create a global variable
cache = Cache()

