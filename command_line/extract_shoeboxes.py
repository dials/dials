#!/usr/bin/env python
#
# dials.extract_shoeboxes.py
#
#  Copyright (C) 2013 Diamond Light Source
#
#  Author: James Parkhurst
#
#  This code is distributed under the BSD license, a copy of which is
#  included in the root directory of this package.

# LIBTBX_SET_DISPATCHER_NAME dev.dials.extract_shoeboxes

from __future__ import division

help_message = '''

This program takes an experiment list and reflections with shoeboxes and
extracts pixels from the images. The shoeboxes are then saved to file.

Examples::

dials.extract_shoeboxes experiments.json reflections.pickle

'''

from libtbx.phil import parse
phil_scope = parse('''
  padding = 0
    .type = int(value_min=0)
    .help = "Add padding around shoebox"

  output {
    reflections = 'shoeboxes.pickle'
      .type = str
      .help = "The integrated output filename"
  }
''')

class Script(object):
  ''' Class to run the script. '''

  def __init__(self):
    ''' Initialise the script. '''

    from dials.util.options import OptionParser
    import libtbx.load_env

    # The script usage
    usage  = "usage: %s [options] experiment.json" % libtbx.env.dispatcher_name

    # Create the parser
    self.parser = OptionParser(
      usage=usage,
      phil=phil_scope,
      epilog=help_message,
      read_experiments=True,
      read_reflections=True,
      read_datablocks=True)

  def run(self):
    ''' Extract the shoeboxes. '''
    from dials.util.options import flatten_reflections
    from dials.util.options import flatten_experiments
    from dials.util.options import flatten_datablocks
    from dials.util import log
    from dials.array_family import flex
    from libtbx.utils import Sorry
    from logging import info

    # Parse the command line
    params, options = self.parser.parse_args(show_diff_phil=False)

    # Configure logging
    log.config()

    # Log the diff phil
    diff_phil = self.parser.diff_phil.as_str()
    if diff_phil is not '':
      info('The following parameters have been modified:\n')
      info(diff_phil)

    # Get the data
    reflections = flatten_reflections(params.input.reflections)
    experiments = flatten_experiments(params.input.experiments)
    datablocks = flatten_datablocks(params.input.datablock)
    if (len(experiments) != 0 and len(datablocks) != 0):
      raise Sorry('Both experiment list and datablocks set')
    elif len(experiments) > 1:
      raise Sorry('More than 1 experiment set')
    elif len(datablocks) > 1:
      raise Sorry('More than 1 datablock set')
    elif len(experiments) == 1:
      imageset = experiments[0].imageset
    elif len(datablocks) == 1:
      imagesets = datablocks[0].extract_imagesets()
      if len(imagesets) != 1:
        raise Sorry('Need 1 imageset, got %d' % len(imagesets))
      imageset = imagesets[0]
    if len(reflections) != 1:
      raise Sorry('Need 1 reflection table, got %d' % len(reflections))
    else:
      reflections = reflections[0]

    # Check the reflections contain the necessary stuff
    assert("bbox" in reflections)
    assert("panel" in reflections)

    # Get some models
    detector = imageset.get_detector()
    scan = imageset.get_scan()
    frame0, frame1 = scan.get_array_range()

    # Add some padding but limit to image volume
    if params.padding > 0:
      info('Adding %d pixels as padding' % params.padding)
      x0, x1, y0, y1, z0, z1 = reflections['bbox'].parts()
      x0 -= params.padding
      x1 += params.padding
      y0 -= params.padding
      y1 += params.padding
      z0 -= params.padding
      z1 += params.padding
      panel = reflections['panel']
      for i in range(len(reflections)):
        width, height = detector[panel[i]].get_image_size()
        if z0[i] < frame0: z0[i] = frame0
        if z1[i] > frame1: z1[i] = frame1
      reflections['bbox'] = flex.int6(x0, x1, y0, y1, z0, z1)

    # Save the old shoeboxes
    if "shoebox" in reflections:
      old_shoebox = reflections['shoebox']
    else:
      old_shoebox = None

    # Allocate the shoeboxes
    reflections["shoebox"] = flex.shoebox(
      reflections["panel"],
      reflections["bbox"],
      allocate=True)

    # Extract the shoeboxes
    reflections.extract_shoeboxes(imageset, verbose=True)

    # Preserve masking
    if old_shoebox is not None:
      info("Applying old shoebox mask")
      new_shoebox = reflections['shoebox']
      for i in range(len(reflections)):
        bbox0 = old_shoebox[i].bbox
        bbox1 = new_shoebox[i].bbox
        mask0 = old_shoebox[i].mask
        mask1 = new_shoebox[i].mask
        mask2 = flex.int(mask1.accessor(), 0)
        x0 = bbox0[0] - bbox1[0]
        x1 = bbox0[1] - bbox0[0] + x0
        y0 = bbox0[2] - bbox1[2]
        y1 = bbox0[3] - bbox0[2] + y0
        z0 = bbox0[4] - bbox1[4]
        z1 = bbox0[5] - bbox0[4] + z0
        mask2[z0:z1,y0:y1,x0:x1] = mask0
        mask1 = mask1.as_1d() | mask2.as_1d()
        mask1.reshape(new_shoebox[i].mask.accessor())
        new_shoebox[i].mask = mask1

    # Saving the reflections to disk
    filename = params.output.reflections
    info('Saving %d reflections to %s' % (len(reflections), filename))
    reflections.as_pickle(filename)

if __name__ == '__main__':
  from dials.util import halraiser
  try:
    script = Script()
    script.run()
  except Exception as e:
    halraiser(e)
