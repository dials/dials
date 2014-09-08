#
# spot_finder.py
#
#  Copyright (C) 2013 Diamond Light Source
#
#  Author: James Parkhurst
#
#  This code is distributed under the BSD license, a copy of which is
#  included in the root directory of this package.

from __future__ import division
from dials.algorithms.peak_finding.threshold import XDSThresholdStrategy


class Extract(object):
  ''' Class to extract a batch of images '''

  def __init__(self, imageset, threshold_image, mask):
      ''' Initialise with imageset and threshold function, both need to be
      picklable for this to be called using multiprocessing. '''
      self.threshold_image = threshold_image
      self.imageset = imageset
      self.mask = mask
      if self.mask is not None:
        detector = self.imageset.get_detector()
        assert(len(self.mask) == len(detector))

  def __call__(self, index):
      ''' Extract pixels from a block of images. '''
      from dials.model.data import PixelList
      from dxtbx.imageset import ImageSweep

      # Get the starting z
      if isinstance(self.imageset, ImageSweep):
        startz = self.imageset.get_array_range()[0] + index[0]
      else:
        ind = self.imageset.indices()
        if len(ind) > 1:
           assert(all(i1+1 == i2 for i1, i2 in zip(ind[0:-1], ind[1:-1])))
        startz = ind[index[0]]

      # Create the list of pixel lists
      plists = [PixelList(p.get_image_size()[::-1], startz)
        for p in self.imageset.get_detector()]

      # Get the trusted ranges
      trange = [p.get_trusted_range() for p in self.imageset.get_detector()]

      # Iterate through the range of images
      for image in self.imageset[index[0]:index[1]]:

        # Ensure image is a tuple of images (for multi-panel support)
        if not isinstance(image, tuple):
          image = (image,)

        # Set the mask
        if self.mask is None:
          mask = []
          for tr, im in zip(trange, image):
            mask.append(im > int(tr[0]))
        else:
          assert(len(self.mask) == len(image))
          mask = self.mask

        # Add the images to the pixel lists
        for pl, im, mk in zip(plists, image, mask):
          pl.add_image(im, self.threshold_image.compute_threshold(im, mk))

      # Return the pixel lists
      return plists


class ProgressUpdater(object):
  ''' Class to update the progress bar from multi processing callback. '''

  def __init__(self, total):
    ''' Initialise the progress bar. '''
    from dials.util.command_line import ProgressBar
    self.total = total
    self.num = 0
    self.progress = ProgressBar(title='Extracting strong pixels from images')

  def __call__(self, result):
    ''' Update the progress bar. '''
    self.num += 1
    percent = 100.0 * (self.num+1) / self.total
    self.progress.update(percent)

  def finished(self):
     ''' Finish the progress bar. '''
     self.progress.finished('Extracted strong pixels from images')


class ExtractSpots(object):
  ''' Class to find spots in an image and extract them into shoeboxes. '''

  def __init__(self, threshold_image, mask=None):
    ''' Initialise the class with the strategy

    Params:
        threshold_image The image thresholding strategy

    '''
    # Set the required strategies
    self.threshold_image = threshold_image
    self.mask = mask

  def __call__(self, imageset):
    ''' Find the spots in the imageset

    Params:
        imageset The imageset to process

    Returns:
        The list of spot shoeboxes

    '''
    from dials.util.command_line import Command
    from dials.model.data import PixelList
    from dials.array_family import flex
    from dxtbx.imageset import ImageSweep
    from libtbx import easy_mp
    from dials.util import mp

    # Change the number of processors if necessary
    nproc = mp.nproc
    if nproc > len(imageset):
      nproc = len(imageset)

    # Extract the pixels in blocks of images in parallel
    Command.start("Extracting strong pixels from images (may take a while)")
    pl = easy_mp.parallel_map(
      func=Extract(imageset, self.threshold_image, self.mask),
      iterable=self._calculate_blocks(imageset, nproc),
      processes=nproc,
      method=mp.method,
      preserve_order=True,
      asynchronous=False)
    Command.end("Extracted strong pixels from images")

    # Merge pixel lists into a single list for each panel
    len_pl = sum(len(p) for p in pl)
    Command.start('Merging {0} pixel lists'.format(len_pl))
    pl = [flex.pixel_list(p).merge() for p in zip(*pl)]
    np = sum([len(p.values()) for p in pl])
    Command.end('Merged {0} pixel lists with {1} pixels'.format(len_pl, np))

    # Extract the pixel lists into a list of reflections
    Command.start('Extracting spots')
    shoeboxes = flex.shoebox()
    if isinstance(imageset, ImageSweep):
      twod = False
    else:
      twod = True
    for i, p in enumerate(pl):
      if p.num_pixels() > 0:
        shoeboxes.extend(flex.shoebox(p, i, 0, twod))
    Command.end('Extracted {0} spots'.format(len(shoeboxes)))

    # Return the shoeboxes
    return shoeboxes

  def _calculate_blocks(self, imageset, nblocks):
    ''' Calculate the blocks. '''
    from math import ceil
    blocks = [0]
    imageset_length = len(imageset)
    assert(nblocks <= imageset_length)
    block_length = int(ceil(imageset_length / nblocks))
    for i in range(nblocks):
      frame = (i + 1) * block_length
      if frame > imageset_length:
        frame = imageset_length
      blocks.append(frame)
      if frame == imageset_length:
        break
    assert(all(b > a for a, b in zip(blocks, blocks[1:])))
    return [(i, j) for i, j in zip(blocks[0:-1], blocks[1:])]


class SpotFinder(object):
  ''' A class to do spot finding and filtering. '''

  def __init__(self, find_spots=None, filter_spots=None, scan_range=None):
    ''' Initialise the class. '''

    # Set the spot finding and filter functions
    assert(find_spots != None and filter_spots != None)
    self.find_spots = find_spots
    self.filter_spots = filter_spots
    self.scan_range = scan_range

  def __call__(self, datablock):
    ''' Do the spot finding. '''
    from dials.array_family import flex

    # Loop through all the imagesets and find the strong spots
    reflections = flex.reflection_table()
    for i, imageset in enumerate(datablock.extract_imagesets()):

      # Find the strong spots in the sweep
      print '-' * 80
      print 'Finding strong spots in imageset %d' % i
      print '-' * 80
      table = self._find_in_imageset(imageset)
      table['id'] = flex.size_t(table.nrows(), i)
      reflections.extend(table)

    # Return the reflections
    return reflections

  def _find_in_imageset(self, imageset):
    ''' Do the spot finding. '''
    from dials.array_family import flex
    from dials.util.command_line import Command
    from dxtbx.imageset import ImageSweep

    # Get the max scan range
    if isinstance(imageset, ImageSweep):
      max_scan_range = imageset.get_array_range()
    else:
      max_scan_range = (0, len(imageset))

    # Get list of scan ranges
    if not self.scan_range or self.scan_range[0] is None:
      scan_range = [max_scan_range]
    else:
      scan_range = self.scan_range

    # Get spots from bits of scan
    spots_all = []
    for scan in scan_range:
      j0, j1 = scan
      assert(j1 > j0 and j0 >= max_scan_range[0] and j1 <= max_scan_range[1])
      print '\nFinding spots in image {0} to {1}...'.format(j0, j1)
      if isinstance(imageset, ImageSweep):
        j0 -= imageset.get_array_range()[0]
        j1 -= imageset.get_array_range()[0]
      spots_all.extend(self.find_spots(imageset[j0:j1]))

    # Get the list of shoeboxes
    shoeboxes = flex.shoebox(spots_all)

    # Calculate the spot centroids
    Command.start('Calculating {0} spot centroids'.format(len(shoeboxes)))
    centroid = shoeboxes.centroid_valid()
    Command.end('Calculated {0} spot centroids'.format(len(shoeboxes)))

    # Calculate the spot intensities
    Command.start('Calculating {0} spot intensities'.format(len(shoeboxes)))
    intensity = shoeboxes.summed_intensity()
    Command.end('Calculated {0} spot intensities'.format(len(shoeboxes)))

    # Create the observations
    observed = flex.observation(shoeboxes.panels(), centroid, intensity)

    # Filter the reflections and select only the desired spots
    flags = self.filter_spots(None,
        sweep=imageset,
        observations=observed,
        shoeboxes=shoeboxes)
    observed = observed.select(flags)
    shoeboxes = shoeboxes.select(flags)

    # Return as a reflection list
    return flex.reflection_table(observed, shoeboxes)
