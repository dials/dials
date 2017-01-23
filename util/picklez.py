#!/usr/bin/env python
#
# picklez.py
#
#  Copyright (C) 2013 Diamond Light Source
#
#  Author: Graeme Winter
#
#  This code is distributed under the BSD license, a copy of which is
#  included in the root directory of this package.
#

# picklez format: pklz - a file containing pickled python objects as strings
# major use case: list of many (objects) which we want to store; to be pushed
#                 into one file as pickled strings of chunks of the list.

from __future__ import absolute_import, division

def chunkify_list(lst, chunk_size):
  '''Return as an iterator chunks of the list, of length chunk_size (as
  number of items.)'''

  start = 0

  while start < len(lst):
    yield lst[start:start + chunk_size]
    start += chunk_size

class Picklez:
  '''Class to pickle objects to a zip file or recover them.'''

  def __init__(self, zipfile_name, zipfile_mode):
    import zipfile
    self._zipfile = zipfile.ZipFile(zipfile_name, zipfile_mode)
    return

  def _list(self):
    return self._zipfile.namelist()

  def _put(self, data, dataname):
    self._zipfile.writestr(dataname, data)
    return

  def _get(self, dataname):
    return self._zipfile.open(dataname, 'r').read()

  def _store(self, obj, name):
    import cPickle as pickle
    self._put(pickle.dumps(obj), name)
    return

  def _recover(self, name):
    import cPickle as pickle
    return pickle.loads(self._get(name))

  def store_list_as_chunks(self, lst, chunk_size):

    for j, chnk in enumerate(chunkify_list(lst, chunk_size)):
      chnk_name = 'chnk%08d.pkl' % j
      self._store(chnk, chnk_name)

    return

  def recover_list_from_chunks(self):

    result = []

    for f in sorted(self._list()):
      if not f.startswith('chnk'):
        continue
      if not f.endswith('.pkl'):
        continue
      result += self._recover(f)

    return result

def dump(lst, zipfile_name, chunk_size):
  picklez = Picklez(zipfile_name, 'w')
  picklez.store_list_as_chunks(lst, chunk_size)
  return

def load(zipfile_name):
  picklez = Picklez(zipfile_name, 'r')
  return picklez.recover_list_from_chunks()

def test():
  import string
  import tempfile

  lst = [string.uppercase for j in range(1000)]

  zipfile_name = tempfile.mktemp(suffix = '.pklz')

  dump(lst, zipfile_name, 100)

  lst2 = load(zipfile_name)

  for l, l2 in zip(lst, lst2):
    assert(l == l2)

  import os
  os.remove(zipfile_name)

  print 'OK'

if __name__ == '__main__':
  test()
