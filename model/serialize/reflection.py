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


class Writer(object):
  ''' Wrapper object to write reflection shoeboxes. '''

  def __init__(self, filename):
    ''' Initialise the writer. '''
    from dials.model.serialize import shoebox
    self._writer = shoebox.Writer(filename)

  def close(self):
    ''' Close the writer. '''
    self._writer.close()

  def write(self, experiment, predictions):
    ''' Write the shoeboxes to disk. '''
    self._writer._write_file(predictions, "predictions")
    self._writer.write(experiment.imageset, predictions)


class Reader(object):
  ''' Read reflections from the given file. '''

  def __init__(self, filename, blocks=None):
    ''' Initialise the writer. '''
    from dials.model.serialize import BlockListIndex
    from dials.model.serialize import shoebox
    from dials.array_family import flex, shared
    self._reader = shoebox.Reader(filename)
    self.predictions = self._reader._read_file("predictions")
    self.blocks = blocks
    batches = [zr for zr, name in self._reader.iter_records()]
    self.lookup = BlockListIndex(shared.tiny_int_2(batches))

  def __len__(self):
    ''' Get number of blocks. '''
    if self.blocks is None:
      return len(self._reader)
    else:
      return len(self.blocks)

  def __getitem__(self, index):
    ''' Read a block of data. '''
    from dials.array_family import flex
    if self.blocks is None:
      zr, ind, sbox = self._reader[index]
      return self._build_table(ind, sbox)
    else:
      zcoord = self.predictions['xyzcal.px'].parts()[2]
      zrange = self.blocks[index]
      indices = flex.size_t()
      shoeboxes = flex.basic_shoebox()
      for i in self.lookup.blocks_in_range(zrange):
        zr, ind, sbox = self._reader[i]
        selection = self.lookup.indices_in_range(
          zrange, zcoord.select(ind))
        indices.extend(ind.select(selection))
        shoeboxes.extend(sbox.select(selection))
      return self._build_table(indices, shoeboxes)

  def __iter__(self):
    ''' Get the number of blocks. '''
    for i in range(len(self)):
      yield self[i]

  def _build_table(self, indices, shoeboxes):
    ''' Build a reflection table from the input. '''
    from dials.array_family import flex
    reflections = self.predictions.select(indices)
    reflections['shoebox'] = flex.shoebox(shoeboxes)
    return reflections


class Extractor(object):
  ''' Helper class to extract data. '''

  def __init__(self, experiment, predictions, filename=None):
    ''' Extract and save to file. '''

    # Set the filename
    if filename is None:
      filename = 'extracted.tar'
    self.filename = filename

    # Write the shoeboxes to file
    writer = Writer(self.filename)
    writer.write(experiment, predictions)

  def extract(self, blocks):
    ''' Iterate through blocks and yield. '''

    # Create the reader
    reader = Reader(self.filename, blocks)

    # Yield all the blocks
    for block in reader:
      yield block
