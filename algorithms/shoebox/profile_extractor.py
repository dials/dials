#
# profile_extractor.py
#
#  Copyright (C) 2013 Diamond Light Source
#
#  Author: James Parkhurst
#
#  This code is distributed under the BSD license, a copy of which is
#  included in the root directory of this package.

from __future__ import division


class PartialProfileExtractor(object):
  ''' A class to extract the profiles from the sweep '''

  def __init__(self, sweep):
    ''' Initialise the class with the sweep etc

    Params:
        sweep The sweep to process

    '''
    self.sweep = sweep

  def __call__(self, panels, bboxes):
    ''' Extract the profiles from the sweep

    Params:
        panels The panel numbers
        bboxes The bounding boxes

    Returns:
        The reflection list

    '''
    from dials.util.command_line import ProgressBar
    from dials.algorithms.shoebox import PartialExtractor

    # Set the number of panels
    image = self.sweep[0]
    if isinstance(image, tuple):
      npanels = len(image)
    else:
      npanels = 1

    # Get the z range
    zrange = self.sweep.get_array_range()

    # Create the class to set all the shoebox pixels
    extractor = PartialExtractor(panels, bboxes, zrange, npanels)

    # For each image in the sweep, get the reflections predicted to have
    # been recorded on the image and copy the pixels from the image to
    # the reflection profile image. Update a progress bar as we go along
    progress = ProgressBar(title = "Extracting frames %d -> %d" % zrange)
    first_array_index = self.sweep.get_array_range()[0]
    for index, image in enumerate(self.sweep, start=first_array_index):

      # Ensure the image is a tuple
      if not isinstance(image, tuple):
        image = (image,)

      # Loop through all the images and add to the extractor
      for panel, im in enumerate(image):
        extractor.add_image(panel, index, im)

      # Update the progress bar
      progress.update(100*(index-first_array_index+1) / len(self.sweep))

    # Get the shoeboxes from the extractor
    shoeboxes = extractor.shoeboxes()
    indices = extractor.shoebox_indices()

    # Finish the progress bar and return the profiles
    progress.finished("Extracted %d profiles from frames %d -> %d"
        % ((len(shoeboxes),) + zrange))

    # Return the indices and shoeboxes
    return (indices, shoeboxes)


class ProfileBlockExtractor(object):
  ''' A class to extract reflections and get them in blocks. '''

  def __init__(self, sweep, predicted, nblocks, filename=None,
               mask=None, gain=None, dark=None, reader=None):
    ''' Initialise the extractor.

    Extract all the data and save into an intermediate format.

    '''
    from dials.model.serialize import partial_shoebox
    from dials.array_family import flex
    from dials.util.command_line import Command

    # Get the gain dark and masks
    self.gain, self.dark, self.mask = self._get_gain_dark_and_mask(
        sweep, gain, dark, mask)

    # Calculate the blocks
    self.blocks = self._calculate_blocks(sweep, nblocks)

    # Create the writer
    if reader == None:
      if filename == None:
        filename = 'extracted.tar'
      writer = partial_shoebox.Writer(filename, predicted, self.blocks)

      # Extract the frames and write to an intermediate format
      panels = flex.size_t([p.panel_number for p in predicted])
      bboxes = flex.int6([p.bounding_box for p in predicted])
      for i, (b0, b1) in enumerate(zip(self.blocks[:-1], self.blocks[1:])):
        extract = PartialProfileExtractor(sweep[b0:b1])
        writer.write(i, *extract(panels, bboxes))
      writer.close()

      # Open the reader ready to get the blocks
      self._reader = partial_shoebox.Reader(filename)
    else:
      assert(len(reader) == len(self))
      self._reader = reader

  def __len__(self):
    ''' Get the number of blocks. '''
    return len(self.blocks) - 1

  def __getitem__(self, index):
    ''' Get the reflections for a particular block. '''
    return self.extract(index)

  def __iter__(self):
    ''' Iterate through all the blocks. '''
    for i in range(len(self)):
      yield self.extract(i)

  def predictions(self):
    ''' Return the list of predictions. '''
    return self._reader.predictions()

  def extract(self, index):
    ''' Get the reflections for a particular block. '''
    from dials.util.command_line import Command
    from dials.array_family import flex
    Command.start('Extracting block %d' % index)
    ind, sb = self._reader.read(index)
    sb = flex.shoebox(sb, self.gain, self.dark, self.mask)
    Command.end('Extracted %d profiles from block %d' % (len(ind), index))
    return ind, sb

  def zrange(self, index):
    ''' Get the frame range of the block '''
    return tuple(self.blocks[index:index+1])

  def _calculate_blocks(self, sweep, nblocks):
    ''' Calculate the blocks. '''
    from math import ceil
    blocks = [0]
    sweep_length = len(sweep)
    assert(nblocks <= sweep_length)
    block_length = int(ceil(sweep_length / nblocks))
    for i in range(nblocks):
      frame = (i + 1) * block_length
      if frame > sweep_length:
        frame = sweep_length
      blocks.append(frame)
      if frame == sweep_length:
        break
    assert(all(b > a for a, b in zip(blocks, blocks[1:])))
    return blocks

  def _get_gain_dark_and_mask(self, sweep, gain, dark, mask):
    ''' Helper function to get gain, dark and mask maps. '''
    from dials.array_family import flex

    # Ensure image is a tuple
    image = sweep[0]
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

class ReflectionBlockExtractor(object):

  def __init__(self, sweep, crystal, predicted, n_sigma, n_blocks,
      filter_by_zeta=0, reader=None, sigma_b=None, sigma_m=None):
    ''' Initialise and extract the reflections. '''
    from dials.algorithms.shoebox import BBoxCalculator
    from dials.util.command_line import Command
    from dials.algorithms import shoebox
    from dials.algorithms import filtering
    from dials.array_family import flex

    if reader == None:

      # These are arrays the length of the sweep
      if sigma_b is not None and sigma_m is not None:
        assert(len(sigma_b) == len(sweep))
        assert(len(sigma_m) == len(sweep))
        compute_bbox = BBoxCalculator(
          sweep.get_beam(), sweep.get_detector(),
          sweep.get_goniometer(), sweep.get_scan(),
          n_sigma * sigma_b, n_sigma * sigma_m)
      else:

        # Create the bbox calculator
        compute_bbox = BBoxCalculator(
          sweep.get_beam(), sweep.get_detector(),
          sweep.get_goniometer(), sweep.get_scan(),
          n_sigma * sweep.get_beam().get_sigma_divergence(deg=False),
          n_sigma * crystal.get_mosaicity(deg=False))

      # Calculate the bounding boxes of all the reflections
      Command.start('Calculating bounding boxes')
      s1 = flex.vec3_double([ r.beam_vector for r in predicted ])
      angle = flex.double([ r.rotation_angle for r in predicted ])
      panel = flex.size_t([ r.panel_number for r in predicted ])
      bbox = compute_bbox(s1, angle, panel)
      for b, r in zip(bbox, predicted):
        r.bounding_box = b
      Command.end('Calculated {0} bounding boxes'.format(len(predicted)))

      # Set all reflections which overlap bad pixels to zero
      Command.start('Filtering reflections by detector mask')
      array_range = sweep.get_scan().get_array_range()
      filtering.by_detector_mask(predicted, sweep[0] >= 0, array_range)
      Command.end('Filtered {0} reflections by detector mask'.format(
          len([r for r in predicted if r.is_valid()])))

      # Filter the reflections by zeta
      if filter_by_zeta > 0:
        Command.start('Filtering reflections by zeta >= {0}'.format(
            filter_by_zeta))
        filtering.by_zeta(sweep.get_goniometer(), sweep.get_beam(),
            predicted, filter_by_zeta)
        Command.end('Filtered {0} reflections by zeta >= {1}'.format(
            len([r for r in predicted if r.is_valid()]), filter_by_zeta))

      # Get only those reflections which are valid
      predicted = predicted.select(predicted.is_valid())

      # Find overlapping reflections
      Command.start('Finding overlapping reflections')
      overlaps = shoebox.find_overlapping(predicted)
      Command.end('Found {0} overlaps'.format(len(overlaps)))

    # Create the profile block extractor
    self._extractor = ProfileBlockExtractor(sweep, predicted, n_blocks, reader=reader)

    # Get the parameters
    delta_d = n_sigma * sweep.get_beam().get_sigma_divergence(deg=False)
    delta_m = n_sigma * crystal.get_mosaicity(deg=False)

    # Create the function to mask the shoebox profiles
    self._mask_profiles = shoebox.Masker(sweep, delta_d, delta_m)

  def __getitem__(self, index):
    ''' Extract a reflection block. '''
    return self.extract(index)

  def __len__(self):
    ''' Get the number of blocks. '''
    return len(self._extractor)

  def __iter__(self):
    ''' Iterate through the blocks. '''
    for i in range(len(self)):
      yield self.extract(i)

  def extract(self, index):
    ''' Extract a block of reflections. '''
    indices, shoeboxes = self._extractor.extract(index)
    reflections = self._extractor.predictions().select(indices)
    for r, s in zip(reflections, shoeboxes):
      r.shoebox = s.data
      r.shoebox_mask = s.mask
      r.shoebox_background = s.background
    self._mask_profiles(reflections, None)
    return reflections
