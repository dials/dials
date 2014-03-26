#!/usr/bin/env python
#
# dials.model.serialize.shoebox.py
#
#  Copyright (C) 2014 Diamond Light Source
#
#  Author: James Parkhurst
#
#  This code is distributed under the BSD license, a copy of which is
#  included in the root directory of this package.

from __future__ import division


#class Reader(object):
  #''' A class to read shoeboxes. '''

  #def __init__(self, filename):
    #''' Load the file and read the index. '''
    #import tarfile
    #import cPickle as pickle
    #self._tar = tarfile.open(filename)
    #index = pickle.load(self._tar.extractfile('index.p'))
    #self._predicted = index['predicted']
    #self._paths = index['paths']
    #self._zrange = index['zrange']

  #def __del__(self):
    #''' close the file. '''
    #self.close()

  #def __getitem__(self, index):
    #''' Read a block of shoeboxes '''
    #return self.read(index)

  #def __len__(self):
    #''' Get the number of blocks. '''
    #return len(self._paths)

  #def __iter__(self):
    #''' Iterate through the blocks. '''
    #for i in range(len(self)):
      #yield self.read(i)

  #def close(self):
    #''' Close the file. '''
    #if hasattr(self, '_tar') and not self._tar.closed:
      #self._tar.close()

  #def read(self, index):
    #''' Read a block of shoeboxes '''
    #return self._load_block(index)

  #def paths(self, index):
    #''' Get the list of paths for a particular block. '''
    #return iter(self._paths[index])

  #def predictions(self):
    #''' Get the predictions. '''
    #return self._predicted

  #def zrange(self):
    #''' Get the z range. '''
    #return self._zrange

  #def _load_block(self, index):
    #''' Load a block of shoeboxes. '''
    #from dials.array_family import flex
    #shoeboxes = flex.partial_shoebox()
    #indices = flex.size_t()
    #for p in self.paths(index):
      #i, s = self._load_partial(p)
      #indices.extend(i)
      #shoeboxes.extend(s)
    #return shoeboxes.merge_all(indices, self._zrange)

  #def _load_partial(self, path):
    #''' Load a list of partial shoeboxes '''
    #import cPickle as pickle
    #return pickle.load(self._tar.extractfile(path))


class Writer(object):
  ''' A class to write shoeboxes '''

  def __init__(self, filename):
    ''' Open the file for writing. '''
    from collections import OrderedDict
    import tarfile
    self._tar = tarfile.open(filename, 'w')

  def __del__(self):
    ''' Close the file. '''
    self.close()

  def close(self):
    ''' Dump the index and close the file. '''
    if hasattr(self, '_tar') and not self._tar.closed:
      self._tar.close()

  def write(self, imageset, predictions):
    ''' Write the shoeboxes for all the predictions. '''
    from dials.util.command_line import ProgressBar
    from dials.model.serialize import BlockList
    zrange = imageset.get_array_range()
    blocks = BlockList(predictions['panel'], predictions['bbox'], zrange)
    p = ProgressBar('Extracting shoeboxes')
    for z, image in enumerate(imageset, start=zrange[0]):
      block = blocks.next(image)
      self._write_block(*block)
      p.update(100.0*(z - zrange[0])/len(imageset))
    p.finished('Extracted shoeboxes')
    assert(blocks.empty())

  def _write_block(self, zrange, indices, shoeboxes):
    ''' Write a block of shoeboxes. '''
    filename = '[%d,%d]' % zrange
    self._write_file((zrange, indices, shoeboxes), filename)

  def _write_file(self, obj, filename):
    ''' Dump a list of partial shoeboxes to file. '''
    import cPickle as pickle
    from StringIO import StringIO
    from time import time
    data = StringIO(pickle.dumps(obj, protocol=pickle.HIGHEST_PROTOCOL))
    info = self._tar.tarinfo()
    info.name = filename
    info.size = data.len
    info.mtime = int(time())
    self._tar.addfile(info, data)
