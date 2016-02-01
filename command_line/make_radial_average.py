#!/usr/bin/env python
#
# dials.make_radial_average.py
#
#  Copyright (C) 2013 Diamond Light Source
#
#  Author: James Parkhurst
#
#  This code is distributed under the BSD license, a copy of which is
#  included in the root directory of this package.

# LIBTBX_SET_DISPATCHER_NAME dev.dials.make_radial_average

from __future__ import division

help_message = '''

This program averages images and makes a radial average over resolution shells

Examples::

dev.dials.make_radial_average experiments.json

'''

# Create the phil scope
from libtbx.phil import parse
phil_scope = parse(
'''

  output {
    filename = 'table.txt'
      .type = str
      .help = "The filename for the output table"
  }

  scan_range = None
    .type = ints(size=2)
    .help = "The scan range to do the average over"

  num_bins = None
    .type = int(value_min=1)
    .help = "The number of bins (default = w+h of image)"

  d_min = None
    .type = float
    .help = "The high resolution limit"

  d_max = None
    .type = float
    .help = "The low resolution limit"

''')


class Script(object):
  ''' The integration program. '''

  def __init__(self):
    '''Initialise the script.'''
    from dials.util.options import OptionParser
    import libtbx.load_env

    # The script usage
    usage  = "usage: %s [options] experiment.json" % libtbx.env.dispatcher_name

    # Create the parser
    self.parser = OptionParser(
      usage=usage,
      phil=phil_scope,
      epilog=help_message,
      read_experiments=True)

  def run(self):
    ''' Perform the integration. '''
    from dials.util.command_line import heading
    from dials.util.options import flatten_reflections, flatten_experiments
    from dials.util import log
    from logging import info, debug
    from time import time
    from libtbx.utils import Sorry
    from dials.array_family import flex

    # Parse the command line
    params, options = self.parser.parse_args(show_diff_phil=False)
    experiments = flatten_experiments(params.input.experiments)
    if len(experiments) != 1:
      self.parser.print_help()
      return

    # Configure logging
    log.config()

    # The imageset
    imageset = experiments[0].imageset

    # Set the scan range
    if params.scan_range is None:
      scan_range = (0, len(imageset))
    else:
      scan_range = params.scan_range
      i0, i1 = scan_range
      if i0 < 0 or i1 > len(imageset):
        raise RuntimeError('Scan range outside image range')
      if i0 >= i1:
        raise RuntimeError('Invalid scan range')

    summed_data = None
    summed_mask = None

    # Loop through images
    for i in range(*scan_range):
      info("Reading image %d" % i)

      # Read image
      data = imageset.get_raw_data(i)
      mask = imageset.get_mask(i)
      assert isinstance(data, tuple)
      assert isinstance(mask, tuple)

      if summed_data is None:
        summed_mask = mask
        summed_data = data
      else:
        summed_data = [ sd + d for sd, d in zip(summed_data, data) ]
        summed_mask = [ sm & m for sm, m in zip(summed_mask, mask) ]

    # Compute min and max and num
    if params.num_bins is None:
      detector = experiments[0].detector
      num_bins = sum(sum(p.get_image_size()) for p in detector)
    if params.d_max is None:
      vmin = 0
    else:
      vmin = (1.0 / d_max)**2
    if params.d_min is None:
      beam = experiments[0].beam
      detector = experiments[0].detector
      params.d_min = detector.get_max_resolution(beam.get_s0())
    vmax = (1.0 / params.d_min)**2

    # Print some info
    info("Min 1/d^2: %f" % vmin)
    info("Max 1/d^2: %f" % vmax)
    info("Num bins:  %d" % num_bins)

    # Compute the radial average
    from dials.algorithms.background import RadialAverage
    radial_average = RadialAverage(beam, detector, vmin, vmax, num_bins)
    for d, m in zip(summed_data, summed_mask):
      radial_average.add(d.as_double(), m)
    mean = radial_average.mean()
    reso = radial_average.inv_d2()

    info("Writing to %s" % params.output.filename)
    with open(params.output.filename, "w") as outfile:
      for r, m in zip(reso, mean):
        outfile.write("%f, %f\n" % (r, m))

if __name__ == '__main__':
  from dials.util import halraiser
  try:
    script = Script()
    script.run()
  except Exception as e:
    halraiser(e)
