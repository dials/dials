#!/usr/bin/env python
#
# dials.model.serialize.partial_shoebox.py
#
#  Copyright (C) 2013 Diamond Light Source
#
#  Author: James Parkhurst
#
#  This code is distributed under the BSD license, a copy of which is
#  included in the root directory of this package.

from __future__ import division


class Reader(object):
    ''' A class to read partial shoeboxes. '''

    def __init__(self, filename):
        ''' Load the file and read the index. '''
        import tarfile
        import cPickle as pickle
        self._tar = tarfile.open(filename)
        index = pickle.load(self._tar.extractfile('index.p'))
        self._predicted = index['predicted']
        self._paths = index['paths']
        self._zrange = index['zrange']

    def __del__(self):
        ''' close the file. '''
        self.close()

    def __getitem__(self, index):
        ''' Read a block of shoeboxes '''
        return self.read(index)

    def __len__(self):
        ''' Get the number of blocks. '''
        return len(self._paths)

    def __iter__(self):
        ''' Iterate through the blocks. '''
        for i in range(len(self)):
            yield self.read(i)

    def close(self):
        ''' Close the file. '''
        if hasattr(self, '_tar') and not self._tar.closed:
            self._tar.close()

    def read(self, index):
        ''' Read a block of shoeboxes '''
        return self._load_block(index)

    def paths(self, index):
        ''' Get the list of paths for a particular block. '''
        return iter(self._paths[index])

    def predictions(self):
        ''' Get the predictions. '''
        return self._predicted

    def zrange(self):
        ''' Get the z range. '''
        return self._zrange

    def _load_block(self, index):
        ''' Load a block of shoeboxes. '''
        from dials.array_family import flex
        shoeboxes = flex.partial_shoebox()
        indices = flex.size_t()
        for p in self.paths(index):
            i, s = self._load_partial(p)
            indices.extend(i)
            shoeboxes.extend(s)
        return shoeboxes.merge_all(indices, self._zrange)

    def _load_partial(self, path):
        ''' Load a list of partial shoeboxes '''
        import cPickle as pickle
        return pickle.load(self._tar.extractfile(path))


class Writer(object):
    ''' A class to write partial shoeboxes '''

    def __init__(self, filename, predicted, blocks):
        ''' Open the file for writing. '''
        from collections import defaultdict
        import tarfile

        # Check the input
        assert(len(blocks) > 1)
        assert(blocks[0] == 0)
        assert(all(b > a for a, b in zip(blocks, blocks[1:])))

        # Set the data members
        self._predicted = predicted
        self._blocks = blocks
        self._tar = tarfile.open(filename, 'w')
        self._paths = defaultdict(list)

        # Create the frame lookup table
        self._lookup = []
        count = 0
        for i in range(max(blocks)):
            if i >= blocks[count+1]:
                count += 1
            self._lookup.append(count)

    def __del__(self):
        ''' Close the file. '''
        self.close()

    def __setitem__(self, index, item):
        ''' Write shoeboxes for a block. '''
        self.write(index, *item)

    def __len__(self):
        ''' Get the number of blocks. '''
        return len(self._blocks) - 1

    def close(self):
        ''' Dump the index and close the file. '''
        if hasattr(self, '_tar') and not self._tar.closed:
            self._dump_index()
            self._tar.close()

    def write(self, index, indices, shoeboxes):
        ''' Write a block of shoeboxes. '''
        self._dump_block(index, indices, shoeboxes)

    def block(self, frame):
        ''' Get the block at a certain frame. '''
        return self._lookup[frame]

    def _dump_block(self, index, indices, shoeboxes):
        ''' Dump a block of shoeboxes to file. '''
        from scitbx.array_family import flex
        assert(len(indices) == len(shoeboxes))
        for block, sind in self._split(indices).iteritems():
            sind = flex.size_t(sind)
            self._dump_shoeboxes(block, indices.select(sind),
                shoeboxes.select(sind))

    def _split(self, indices):
        ''' Split the shoeboxes based on the target block. '''
        from collections import defaultdict
        sperb = defaultdict(list)
        for i in range(len(indices)):
            frame = int(self._predicted[indices[i]].frame_number)
            sperb[self.block(frame)].append(i)
        return sperb

    def _dump_shoeboxes(self, block, indices, shoeboxes):
        ''' Dump a list of partial shoeboxes to file. '''
        import cPickle as pickle
        from StringIO import StringIO
        from uuid import uuid4
        from time import time
        data = StringIO(pickle.dumps((indices, shoeboxes),
            protocol=pickle.HIGHEST_PROTOCOL))
        info = self._tar.tarinfo()
        info.name = '%s.p' % uuid4().hex
        info.size = data.len
        info.mtime = int(time())
        self._tar.addfile(info, data)
        self._paths[block].append(info.name)

    def _dump_index(self):
        ''' Dump the index to the file. '''
        import cPickle as pickle
        from StringIO import StringIO
        from time import time
        index = {
          'paths' : self._paths,
          'predicted' : self._predicted,
          'zrange' : (min(self._blocks), max(self._blocks))
        }
        data = StringIO(pickle.dumps(index, protocol=pickle.HIGHEST_PROTOCOL))
        info = self._tar.tarinfo()
        info.name = 'index.p'
        info.size = data.len
        info.mtime = int(time())
        self._tar.addfile(info, data)
