#!/usr/bin/env python
#
#  reflection_block.py
#
#  Copyright (C) 2013 Diamond Light Source
#
#  Author: James Parkhurst
#
#  This code is distributed under the BSD license, a copy of which is
#  included in the root directory of this package.

from __future__ import division


class ReflectionBlockExtractor(object):
  ''' A class to extract blocks of reflections. '''

  def __init__(self, filename, block_size,
               imageset, reflections=None,
               gain=None, dark=None, mask=None):
    ''' Initialise the extractor. '''
    from dials.model.serialize import extract_shoeboxes_to_file
    from dials.model.serialize import ShoeboxBlockImporter
    from dials.array_family import flex
    import cPickle as pickle
    from math import log10, floor

    # Save the imageset
    self.imageset = imageset

    # Ensure we have lookup tables
    gain, dark, mask = self._check_lookup(imageset, gain, dark, mask)

    # Reorder the reflections and extract shoeboxes
    if reflections:
      self.reflections = self._reorder_reflections(reflections)
      extract_shoeboxes_to_file(filename, imageset, reflections)

    # Calculate the blocks in images
    if imageset.get_scan() is None:
      self._blocks = [0,len(imageset)]
    else:
      self._blocks = self._compute_blocks(imageset.get_scan(), block_size)
    print "Extracting reflections from the following blocks of images:"
    npad = int(floor(log10(max(self._blocks)))) + 1
    format_string = ' %%%dd -> %%%dd' % (npad, npad)
    for i in range(len(self._blocks)-1):
      print format_string % (self._blocks[i], self._blocks[i+1])

    # Construct the importer
    if gain and dark and mask:
      self._importer = ShoeboxBlockImporter(
        filename, flex.size_t(self._blocks), gain, dark, mask)
    else:
      self._importer = ShoeboxBlockImporter(
        filename, flex.size_t(self._blocks))

    # If reflections not set, then read from blob
    if reflections is None:
      self.reflections = pickle.loads(self._importer.blob())

  def block(self, index):
    ''' Get the block. '''
    return tuple(self._blocks[index:index+1])

  def __len__(self):
    ''' Return the number of blocks. '''
    return len(self._importer)

  def __getitem__(self, index):
    ''' Get a block of reflections. '''
    from dials.util.command_line import Command

    # Get the indices and shoeboxes
    Command.start('Extracting block %d' % index)
    indices, shoeboxes = self._importer[index]

    # Create the partial array of reflections
    reflections = self.reflections.select(indices)
    reflections['shoebox'] = shoeboxes
    nref = len(reflections)
    Command.end('Extracted %d reflections from block %d' % (nref, index))

    # return the indices and reflections
    return indices, reflections

  def __iter__(self):
    ''' Iterate through the blocks. '''
    for i in range(len(self)):
      yield self[i]

  def _reorder_reflections(self, reflections):
    ''' Filter the reflections and sort them by z. '''
    from dials.array_family import flex

    # Sort the reflections by z
    z = reflections['xyzcal.px'].parts()[2]
    indices = flex.size_t(sorted(range(len(z)), key=lambda x: z[x]))
    reflections.reorder(indices)

    # Return the reflections
    return reflections

  def _compute_blocks(self, scan, block_size):
    ''' Compute the number of blocks. '''
    from math import ceil
    phi0, dphi = scan.get_oscillation(deg=True)
    nframes = scan.get_num_images()
    frame0 = scan.get_array_range()[0]
    assert(block_size >= dphi)
    block_length = float(block_size) / dphi
    nblocks = int(ceil(nframes / block_length))
    assert(nblocks <= nframes)
    block_length = int(ceil(nframes / nblocks))
    blocks = [frame0]
    for i in range(nblocks):
      frame = frame0 + (i + 1) * block_length
      if frame > frame0+nframes:
        frame = frame0+nframes
      blocks.append(frame)
      if frame == frame0+nframes:
        break
    assert(all(b > a for a, b in zip(blocks, blocks[1:])))
    return blocks

  def _check_lookup(self, imageset, gain, dark, mask):
    ''' Ensure the lookup tables are Ok '''
    from dials.array_family import flex
    detector = imageset.get_detector()
    npanels = len(detector)
    if gain is None:
      gain = []
      for i in range(npanels):
        isize = detector[i].get_image_size()[::-1]
        gain.append(flex.double(flex.grid(isize),1.0))
    if dark is None:
      dark = []
      for i in range(npanels):
        isize = detector[i].get_image_size()[::-1]
        dark.append(flex.double(flex.grid(isize),0.0))
    if mask is None:
      mask = []
      image = imageset[0]
      if npanels == 1:
        image = (image,)
      for i in range(npanels):
        trusted_range = detector[i].get_trusted_range()
        m = image[i] > int(trusted_range[0])
        mask.append(m)
    return tuple(gain), tuple(dark), tuple(mask)
