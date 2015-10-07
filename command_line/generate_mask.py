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

  untrusted
    .multiple = True
  {

    panel = 0
      .type = int
      .help = "The panel number"

    circle = None
      .type = ints(3)
      .help = "An untrusted circle (xc, yc, r)"

    rectangle = None
      .type = ints(4)
      .help = "An untrusted rectangle (x0, x1, y0, y1)"
  }

  output {
    mask = mask.pickle
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
    from math import floor, ceil

    # Get the detector
    detector = imageset.get_detector()

    # Get the first image
    image = imageset[0]
    if not isinstance(image, tuple):
      image = (image,)
    assert(len(detector) == len(image))

    # Create the mask for each image
    masks = []
    for index, (im, panel) in enumerate(zip(image, detector)):

      # The image width height
      height, width = im.all()

      # Create the basic mask from the trusted range
      if self.params.use_trusted_range:
        low, high = panel.get_trusted_range()
        imd = im.as_double()
        mask = (imd > low) & (imd < high)
      else:
        mask = flex.bool(flex.grid(im.all()), True)

      # Add a border around the image
      if self.params.border > 0:
        print "Generating border mask:"
        print " border = %d" % self.params.border
        border = self.params.border
        height, width = mask.all()
        borderx = flex.bool(flex.grid(border, width), False)
        bordery = flex.bool(flex.grid(height, border), False)
        mask[0:border,:] = borderx
        mask[-border:,:] = borderx
        mask[:,0:border] = bordery
        mask[:,-border:] = bordery

      # Apply the untrusted region
      for region in self.params.untrusted:
        if region.panel == index:
          if region.circle is not None:
            xc, yc, radius = region.circle
            x0 = int(floor(xc - radius))
            y0 = int(floor(yc - radius))
            x1 = int(ceil(xc + radius))
            y1 = int(ceil(yc + radius))
            assert x1 > x0
            assert y1 > y0
            assert x0 >= 0
            assert y0 >= 0
            assert x1 <= width
            assert y1 <= height
            print "Generating circle mask:"
            print " panel = %d" % region.panel
            print " xc = %d" % xc
            print " yc = %d" % yc
            print " radius = %d" % radius
            xc -= x0
            yc -= y0
            r2 = radius * radius
            circ = flex.bool(flex.grid(y1-y0,x1-x0),False)
            for j in range(y1-y0):
              for i in range(x1-x0):
                if (i - xc)**2 + (j - yc)**2 > r2:
                  circ[j,i] = True
            mask[y0:y1,x0:x1] = mask[y0:y1,x0:x1] & circ
          if region.rectangle is not None:
            x0, x1, y0, y1 = region.rectangle
            assert x1 > x0
            assert y1 > y0
            assert x0 >= 0
            assert y0 >= 0
            assert x1 <= width
            assert y1 <= height
            print "Generating rectangle mask:"
            print " panel = %d" % region.panel
            print " x0 = %d" % x0
            print " y0 = %d" % y0
            print " x1 = %d" % x1
            print " y1 = %d" % y1
            rect = flex.bool(flex.grid(y1-y0,x1-x0),False)
            mask[y0:y1,x0:x1] = rect

      # Add to the list
      masks.append(mask)

    # Return the mask
    return tuple(masks)


class Script(object):
  ''' A class to encapsulate the script. '''

  def __init__(self):
    ''' Initialise the script. '''
    from dials.util.options import OptionParser
    import libtbx.load_env

    # Create the parser
    usage = "usage: %s [options] datablock.json" % libtbx.env.dispatcher_name
    self.parser = OptionParser(
      usage=usage,
      phil=phil_scope,
      read_datablocks=True)

  def run(self):
    ''' Run the script. '''
    from dials.util.options import flatten_datablocks
    from libtbx.utils import Sorry
    import cPickle as pickle

    # Parse the command line arguments
    params, options = self.parser.parse_args(show_diff_phil=True)
    datablocks = flatten_datablocks(params.input.datablock)

    # Check number of args
    if len(datablocks) == 0:
      self.parser.print_help()
      return

    if len(datablocks) != 1:
      raise Sorry('exactly 1 datablock must be specified')
    datablock = datablocks[0]
    imagesets = datablock.extract_imagesets()
    if len(imagesets) != 1:
      raise Sorry('datablock must contain exactly 1 imageset')
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
