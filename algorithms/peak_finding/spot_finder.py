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


class Extract(object):
  '''
  Class to extract a batch of images

  '''

  def __init__(self, imageset, threshold_image, mask):
      '''
      Initialise with imageset and threshold function, both need to be picklable
      for this to be called using multiprocessing.

      :param imageset: The imageset to process
      :param threshold_image: The threshold algorithm
      :param mask: The mask to use

      '''
      self.threshold_image = threshold_image
      self.imageset = imageset
      self.mask = mask
      if self.mask is not None:
        detector = self.imageset.get_detector()
        assert(len(self.mask) == len(detector))

  def __call__(self, index):
      '''
      Extract pixels from a block of images.

      :param index: The image indices to extract
      :return: The extracted spots

      '''
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

      # Iterate through the range of images
      for ind in range(*index):

        # Get the image and mask
        image = self.imageset.get_corrected_data(ind)
        mask = self.imageset.get_mask(ind)

        # Set the mask
        if self.mask is not None:
          assert(len(self.mask) == len(mask))
          mask = tuple(m1 & m2 for m1, m2 in zip(mask, self.mask))

        # Add the images to the pixel lists
        for pl, im, mk in zip(plists, image, mask):
          pl.add_image(im, self.threshold_image.compute_threshold(im, mk))

      # Return the pixel lists
      return plists


class ProgressUpdater(object):
  '''
  Class to update the progress bar from multi processing callback.

  '''

  def __init__(self, total):
    '''
    Initialise the progress bar.

    :param total: The total to count up to

    '''
    from dials.util.command_line import ProgressBar
    self.total = total
    self.num = 0
    self.progress = ProgressBar(title='Extracting strong pixels from images')

  def __call__(self, result):
    '''
    Update the progress bar.

    :param result: The result

    '''
    self.num += 1
    percent = 100.0 * (self.num+1) / self.total
    self.progress.update(percent)

  def finished(self):
     '''
     Finish the progress bar.

     '''
     self.progress.finished('Extracted strong pixels from images')


class ExtractSpots(object):
  '''
  Class to find spots in an image and extract them into shoeboxes.

  '''

  def __init__(self, threshold_image, mask=None,
               mp_method='multiprocessing', nproc=1):
    '''
    Initialise the class with the strategy

    :param threshold_image: The image thresholding strategy
    :param mask: The mask to use
    :param mp_method: The multi processing method
    :param nproc: The number of processors

    '''
    # Set the required strategies
    self.threshold_image = threshold_image
    self.mask = mask
    self.mp_method = mp_method
    self.nproc = nproc

  def __call__(self, imageset):
    '''
    Find the spots in the imageset

    :param imageset: The imageset to process
    :return: The list of spot shoeboxes

    '''
    from dials.util.command_line import Command
    from dials.array_family import flex
    from dxtbx.imageset import ImageSweep
    from libtbx import easy_mp
    from logging import info

    # Change the number of processors if necessary
    nproc = self.nproc
    if nproc > len(imageset):
      nproc = len(imageset)

    # Extract the pixels in blocks of images in parallel
    info("Extracting strong pixels from images (may take a while)")
    pl = easy_mp.parallel_map(
      func=Extract(imageset, self.threshold_image, self.mask),
      iterable=self._calculate_blocks(imageset, nproc),
      processes=nproc,
      method=self.mp_method,
      preserve_order=True,
      asynchronous=False)
    info("Extracted strong pixels from images")

    # Merge pixel lists into a single list for each panel
    len_pl = sum(len(p) for p in pl)
    info('Merging {0} pixel lists'.format(len_pl))
    pl = [flex.pixel_list(p).merge() for p in zip(*pl)]
    np = sum([len(p.values()) for p in pl])
    info('Merged {0} pixel lists with {1} pixels'.format(len_pl, np))

    # Extract the pixel lists into a list of reflections
    info('Extracting spots')
    shoeboxes = flex.shoebox()
    if isinstance(imageset, ImageSweep):
      twod = False
    else:
      twod = True
    for i, p in enumerate(pl):
      if p.num_pixels() > 0:
        shoeboxes.extend(flex.shoebox(p, i, 0, twod))
    info('Extracted {0} spots'.format(len(shoeboxes)))

    # Return the shoeboxes
    return shoeboxes

  def _calculate_blocks(self, imageset, nblocks):
    '''
    Calculate the blocks.

    :param imageset: The imageset
    :param nblocks: The number of blocks to process in
    :return: The block size

    '''
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
  '''
  A class to do spot finding and filtering.

  '''

  def __init__(self, find_spots=None, filter_spots=None, scan_range=None,
               write_hot_mask=True):
    '''
    Initialise the class.

    :param find_spots: The spot finding algorithm
    :param filter_spots: The spot filtering algorithm
    :param scan_range: The scan range to find spots over

    '''

    # Set the spot finding and filter functions
    assert(find_spots != None and filter_spots != None)
    self.find_spots = find_spots
    self.filter_spots = filter_spots
    self.scan_range = scan_range
    self.write_hot_mask = write_hot_mask

  def __call__(self, datablock):
    '''
    Do the spot finding.

    :param datablock: The datablock to process
    :return: The observed spots

    '''
    from dials.array_family import flex
    from logging import info
    import cPickle as pickle

    # Loop through all the imagesets and find the strong spots
    reflections = flex.reflection_table()
    for i, imageset in enumerate(datablock.extract_imagesets()):

      # Find the strong spots in the sweep
      info('-' * 80)
      info('Finding strong spots in imageset %d' % i)
      info('-' * 80)
      table, hot_mask = self._find_in_imageset(imageset)
      table['id'] = flex.size_t(table.nrows(), i)
      reflections.extend(table)
      if imageset.external_lookup.mask.data is not None:
        and_mask = []
        for m1, m2 in zip(imageset.external_lookup.mask.data, hot_mask):
          and_mask.append(m1 & m2)
        imageset.external_lookup.mask.data = tuple(and_mask)
      else:
        imageset.external_lookup.mask.data = hot_mask
      imageset.external_lookup.mask.filename = "hot_mask_%d.pickle" % i

      # Write the hot mask
      if self.write_hot_mask:
        with open(imageset.external_lookup.mask.filename, "wb") as outfile:
          pickle.dump(hot_mask, outfile, protocol=pickle.HIGHEST_PROTOCOL)

    reflections.set_flags(
      flex.size_t_range(len(reflections)), reflections.flags.strong)

    # Return the reflections
    return reflections

  def _find_in_imageset(self, imageset):
    '''
    Do the spot finding.

    :param imageset: The imageset to process
    :return: The observed spots

    '''
    from dials.array_family import flex
    from dials.util.command_line import Command
    from dxtbx.imageset import ImageSweep
    from logging import info

    # Get the max scan range
    if isinstance(imageset, ImageSweep):
      max_scan_range = imageset.get_array_range()
    else:
      max_scan_range = (0, len(imageset))

    # Get list of scan ranges
    if not self.scan_range or self.scan_range[0] is None:
      scan_range = [(max_scan_range[0]+1, max_scan_range[1])]
    else:
      scan_range = self.scan_range

    # Get spots from bits of scan
    spots_all = []
    for scan in scan_range:
      j0, j1 = scan
      assert(j1 >= j0 and j0 > max_scan_range[0] and j1 <= max_scan_range[1])
      info('\nFinding spots in image {0} to {1}...'.format(j0, j1))
      j0 -= 1
      if isinstance(imageset, ImageSweep):
        j0 -= imageset.get_array_range()[0]
        j1 -= imageset.get_array_range()[0]
      spots_all.extend(self.find_spots(imageset[j0:j1]))

    # Get the list of shoeboxes
    shoeboxes = flex.shoebox(spots_all)

    # Calculate the spot centroids
    info('Calculating {0} spot centroids'.format(len(shoeboxes)))
    centroid = shoeboxes.centroid_valid()
    info('Calculated {0} spot centroids'.format(len(shoeboxes)))

    # Calculate the spot intensities
    info('Calculating {0} spot intensities'.format(len(shoeboxes)))
    intensity = shoeboxes.summed_intensity()
    info('Calculated {0} spot intensities'.format(len(shoeboxes)))

    # Create the observations
    observed = flex.observation(shoeboxes.panels(), centroid, intensity)

    # Find spots which cover the whole scan range
    bbox = flex.int6([sbox.bbox for sbox in shoeboxes])
    z0, z1 = bbox.parts()[4:6]
    zr = z1 - z0
    assert zr.all_gt(0)
    possible_hot_spots = (zr == len(imageset))
    num_possible_hot_spots = possible_hot_spots.count(True)
    info('Found %d possible hot spots' % num_possible_hot_spots)

    # Create the hot pixel mask
    hot_mask = tuple(flex.bool(flex.grid(p.get_image_size()[::-1]), True) for p
                     in imageset.get_detector())
    if num_possible_hot_spots > 0:
      hot_shoeboxes = shoeboxes.select(possible_hot_spots)
      for sbox in hot_shoeboxes:
        x0, x1, y0, y1 = sbox.bbox[0:4]
        m = sbox.mask
        p = sbox.panel
        for y in range(m.all()[1]):
          for x in range(m.all()[2]):
            if m[:,y:y+1,x:x+1].all_ne(0):
              hot_mask[p][y0+y,x0+x] = False
    info('Found %d possible hot pixel(s)' % hot_mask.count(False))

    # Filter the reflections and select only the desired spots
    flags = self.filter_spots(None,
        sweep=imageset,
        observations=observed,
        shoeboxes=shoeboxes)
    observed = observed.select(flags)
    shoeboxes = shoeboxes.select(flags)

    # Return as a reflection list
    return flex.reflection_table(observed, shoeboxes), hot_mask
