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
from dials.interfaces.peak_finding import SpotFinderInterface
from dials.algorithms.peak_finding.threshold import XDSThresholdStrategy


class Extract(object):
  ''' Class to extract a batch of images '''

  def __init__(self, sweep, threshold_image):
      ''' Initialise with sweep and threshold function, both need to be
      picklable for this to be called using multiprocessing. '''
      self.threshold_image = threshold_image
      self.sweep = sweep

  def __call__(self, index):
      ''' Extract pixels from a block of images. '''
      from dials.model.data import PixelList

      # Create the list of pixel lists
      plists = [PixelList(p.get_image_size()[::-1], index[0])
        for p in self.sweep.get_detector()]

      # Iterate through the range of images
      for image in self.sweep[index[0]:index[1]]:

        # Ensure image is a tuple of images (for multi-panel support)
        if not isinstance(image, tuple):
          image = (image,)

        # Add the images to the pixel lists
        for pl, im in zip(plists, image):
          pl.add_image(im, self.threshold_image(im))

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

  def __init__(self, threshold_image):
    ''' Initialise the class with the strategy

    Params:
        threshold_image The image thresholding strategy

    '''
    # Set the required strategies
    self.threshold_image = threshold_image

  def __call__(self, sweep):
    ''' Find the spots in the sweep

    Params:
        sweep The sweep to process

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
    if nproc > len(sweep):
      nproc = len(sweep)

    # Extract the pixels in blocks of images in parallel
    progress = ProgressUpdater(nproc)
    pl = easy_mp.parallel_map(
      func=Extract(sweep, self.threshold_image),
      iterable=self._calculate_blocks(sweep, nproc),
      processes=nproc,
      method=mp.method,
      preserve_order=True,
      asynchronous=True,
      callback=progress)
    progress.finished()

    # Merge pixel lists into a single list for each panel
    len_pl = sum(len(p) for p in pl)
    Command.start('Merging {0} pixel lists'.format(len_pl))
    pl = [flex.pixel_list(p).merge() for p in zip(*pl)]
    np = sum([len(p.values()) for p in pl])
    Command.end('Merged {0} pixel lists with {1} pixels'.format(len_pl, np))

    # Extract the pixel lists into a list of reflections
    Command.start('Extracting spots')
    shoeboxes = flex.shoebox()
    if isinstance(sweep, ImageSweep):
      twod = False
    else:
      twod = True
    for i, p in enumerate(pl):
      if p.num_pixels() > 0:
        shoeboxes.extend(flex.shoebox(p, i, 0, twod))
    Command.end('Extracted {0} spots'.format(len(shoeboxes)))

    # Return the shoeboxes
    return shoeboxes

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
    return [(i, j) for i, j in zip(blocks[0:-1], blocks[1:])]


class SpotFinder(object):
  ''' A class to do spot finding and filtering. '''

  def __init__(self, find_spots=None, filter_spots=None, scan_range=None):
    ''' Initialise the class. '''

    # Set the spot finding and filter functions
    assert(find_spots != None and filter_spots != None)
    self.find_spots = find_spots
    self.filter_spots = filter_spots

    # Set the scan range
    self.scan_range = scan_range

  def __call__(self, sweep):
    ''' Do the spot finding '''
    from dials.model.data import ReflectionList
    from dials.array_family import flex
    from dials.util.command_line import Command
    from dxtbx.imageset import ImageSweep

    # Get list of scan ranges
    if not self.scan_range:
      if isinstance(sweep, ImageSweep):
        scan_range = [sweep.get_array_range()]
      else:
        scan_range = [(0, len(sweep))]
    else:
      scan_range = self.scan_range

    # Get spots from bits of scan
    spots_all = []
    for scan in scan_range:
      j0, j1 = scan
      assert(j1 > j0 and j0 >= 0 and j1 <= len(sweep))
      print '\nFinding spots in image {0} to {1}...'.format(j0, j1)
      spots_all.extend(self.find_spots(sweep[j0:j1]))

    # Get the list of shoeboxes
    shoeboxes = flex.shoebox(spots_all)

    # Calculate the spot centroids
    Command.start('Calculating {0} spot centroids'.format(len(shoeboxes)))
    centroid = shoeboxes.centroid_valid()
    Command.end('Calculated {0} spot centroids'.format(len(shoeboxes)))

    # Calculate the spot intensities
    Command.start('Calculating {0} spot intensities'.format(len(shoeboxes)))
    intensity = shoeboxes.summed_intensity_valid()
    Command.end('Calculated {0} spot intensities'.format(len(shoeboxes)))

    # Create the observations
    observed = flex.observation(shoeboxes.panels(), centroid, intensity)

    # Filter the reflections and select only the desired spots
    flags = self.filter_spots(None,
        sweep=sweep,
        observations=observed,
        shoeboxes=shoeboxes)
    observed = observed.select(flags)
    shoeboxes = shoeboxes.select(flags)

    # Return as a reflection list
    return ReflectionList(observed, shoeboxes)
