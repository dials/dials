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

  def __init__(self, filename, blocks=None, mask=None, gain=None, dark=None):
    ''' Initialise the writer. '''
    from dials.model.serialize import BlockListIndex
    from dials.model.serialize import shoebox
    from dials.array_family import flex, shared
    self._reader = shoebox.Reader(filename)
    self.predictions = self._reader._read_file("predictions")
    self.blocks = blocks
    batches = [zr for zr, name in self._reader.iter_records()]
    self.lookup = BlockListIndex(shared.tiny_int_2(batches))
    self.mask = mask
    self.gain = gain
    self.dark = dark

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
    if self.mask is None or self.gain is None or self.dark is None:
      reflections['shoebox'] = flex.shoebox(shoeboxes)
    else:
      reflections['shoebox'] = flex.shoebox(shoeboxes,
        self.mask, self.gain, self.dark)
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

    # Save the experiment
    self.experiment = experiment


  def extract(self, blocks, mask=None, gain=None, dark=None):
    ''' Iterate through blocks and yield. '''

    # Create the reader
    reader = Reader(self.filename, blocks, *self._lookup(mask, gain, dark))

    # Yield all the blocks
    for block in reader:
      yield block

  def _lookup(self, gain, dark, mask):
    ''' Helper function to get gain, dark and mask maps. '''
    from dials.array_family import flex

    # Ensure image is a tuple
    image = self.experiment.imageset[0]
    if not isinstance(image, tuple):
      image = (image,)

    # Get the mask in tuple of masks form
    if mask:
      if not isinstance(mask, tuple):
        mask = (mask,)
    else:
      mask = tuple([im >= 0 for im in image])

    # Get the gain in tuple of gains form
    if gain:
      if not isinstance(gain, tuple):
        gain = (gain,)
    else:
      gain = tuple([flex.double(flex.grid(im.all()), 1) for im in image])

    # Get the dark in tuple of darks form
    if dark:
      if not isinstance(dark, tuple):
        dark = (dark,)
    else:
      dark = tuple([flex.double(flex.grid(im.all()), 0) for im in image])

    return gain, dark, mask
