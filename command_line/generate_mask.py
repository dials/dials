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
from iotbx.phil import parse

phil_scope = parse("""
  border = 0
    .type = int
    .help = "The border around the edge of the image."

  use_trusted_range = True
    .type = bool
    .help = "Use the trusted range to mask bad pixels."

  output {
    mask = mask.p
      .type = str
      .help = "Name of output mask file"
  }
""", process_includes=True)

class MaskGenerator(object):
  ''' Generate a mask. '''

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
      if self.params.use_trusted_range:
        low, high = panel.get_trusted_range()
        imd = im.as_double()
        mask = (imd > low) & (imd < high)
      else:
        mask = flex.bool(flex.grid(im.all()), True)

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


class Script(object):
  ''' A class to encapsulate the script. '''

  def __init__(self):
    ''' Initialise the script. '''
    from dials.util.options import OptionParser

    # Create the parser
    usage = "usage: %prog [options] datablock.json"
    self.parser = OptionParser(
      usage=usage,
      phil=phil_scope,
      read_datablocks=True)

  def run(self):
    ''' Run the script. '''
    from dials.util.options import flatten_datablocks
    import cPickle as pickle

    # Parse the command line arguments
    params, options = self.parser.parse_args(show_diff_phil=True)
    datablocks = flatten_datablocks(params.input.datablock)

    # Check number of args
    if len(datablocks) == 0:
      self.parser.print_help()
      return

    assert(len(datablocks) == 1)
    datablock = datablocks[0]
    imagesets = datablock.extract_imagesets()
    assert(len(imagesets) == 1)
    imageset = imagesets[0]

    # Generate the mask
    generator = MaskGenerator(params)
    mask = generator.generate(imageset)

    # Save the mask to file
    pickle.dump(mask, open(params.output.mask, "w"))

if __name__ == '__main__':
  from dials.util import halraiser
  try:
    script = Script()
    script.run()
  except Exception as e:
    halraiser(e)
