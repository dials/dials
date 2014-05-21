#!/usr/bin/env python
#
# generate_mask.py
#
#  Copyright (C) 2013 Diamond Light Source
#
#  Author: James Parkhurst
#
#  This code is distributed under the BSD license, a copy of which is
#  included in the root directory of this package.

from __future__ import division

class MaskGenerator(object):
  ''' Generate a mask. '''

  phil = '''

    border = 0
      .type = int
      .help = "The border around the edge of the image."

  '''

  def __init__(self, params):
    ''' Set the parameters. '''
    self.params = params

  def generate(self, imageset):
    ''' Generate the mask. '''
    from dials.array_family import flex

    # Get the detector
    detector = imageset.get_detector()

    # Get the first image
    image = imageset[0]
    if not isinstance(image, tuple):
      image = (image,)
    assert(len(detector) == len(image))

    # Create the mask for each image
    masks = []
    for im, panel in zip(image, detector):

      # Create the basic mask from the trusted range
      low, high = panel.get_trusted_range()
      imd = im.as_double()
      mask = (imd > low) & (imd < high)

      # Add a border around the image
      if self.params.border > 0:
        border = self.params.border
        height, width = mask.all()
        borderx = flex.bool(flex.grid(border, width), False)
        bordery = flex.bool(flex.grid(height, border), False)
        mask[0:border,:] = borderx
        mask[-border:,:] = borderx
        mask[:,0:border] = bordery
        mask[:,-border:] = bordery

      # Add to the list
      masks.append(mask)

    # Return the mask
    return tuple(masks)


if __name__ == '__main__':

  from libtbx.phil import parse
  from libtbx.phil import command_line
  from optparse import OptionParser
  from dxtbx.datablock import DataBlockFactory
  import cPickle as pickle

  usage = "usage: %prog [options] datablock.json"
  parser = OptionParser(usage)

  # Parse the command line arguments
  (options, args) = parser.parse_args()

  # Get phil parameters
  master_phil = parse(MaskGenerator.phil)

  # Process the command line arguments
  cmd = command_line.argument_interpreter(master_params = master_phil)
  working_phils = []
  unhandled = []
  for arg in args:
    try:
      working_phils.append(cmd.process_arg(arg))
    except Exception:
      unhandled.append(arg)
  working_phil = master_phil.fetch(sources=working_phils)
  params = working_phil.extract()

  # Check number of args
  if len(unhandled) == 0:
    parser.print_help()
    exit(0)

  # Load the data block
  datablocks = []
  for arg in unhandled:
    datablocks.extend(DataBlockFactory.from_json_file(arg))
  assert(len(datablocks) == 1)
  datablock = datablocks[0]
  imagesets = datablock.extract_imagesets()
  assert(len(imagesets) == 1)
  imageset = imagesets[0]

  # Generate the mask
  generator = MaskGenerator(params)
  mask = generator.generate(imageset)

  # Save the mask to file
  pickle.dump(mask, open("mask.p", "w"))
