#!/usr/bin/env python
#
# perturb_reflection_phi.py
#
#  Copyright (C) (2013) STFC Rutherford Appleton Laboratory, UK.
#
#  Author: David Waterman.
#
#  This code is distributed under the BSD license, a copy of which is
#  included in the root directory of this package.

"""This script applies a sinusoidal perturbation of up to 0.5 degree to the
reflecting phi for a set of input reflections. This is a way to exercise
scan-varying refinement, which should be able to model the sinusoidal
variation."""

from __future__ import division

class ScriptRunner(object):
  """Class to run script."""

  def __init__(self, options, reflections_filename):
    """Setup the script."""

    # Filename data
    self.output_filename = options.output_file
    self.reflections_filename = reflections_filename

  def __call__(self):
    """Run the script."""
    from dials.model.serialize import load, dump
    import cPickle as pickle
    from math import sin, pi

    # Load the reflection list
    print 'Loading reflections from {0}'.format(self.reflections_filename)
    rlist = pickle.load(open(self.reflections_filename, 'r'))

    # Loop through all the reflections
    for r in rlist:

      d_phi = sin(r.rotation_angle) * 0.5 * pi / 180
      r.rotation_angle += d_phi
      r.centroid_position = r.image_coord_mm + (r.rotation_angle, )

    # Write out reflections
    if self.output_filename is not None:

      print 'Saving reflections to {0}'.format(self.output_filename)
      pickle.dump(rlist, open(self.output_filename, 'wb'),
          pickle.HIGHEST_PROTOCOL)

if __name__ == '__main__':

  from optparse import OptionParser

  # Specify the command line options
  usage  = "usage: %prog [options] reflections.pickle"

  # Create an option parser
  parser = OptionParser(usage)

  parser.add_option('-o', '--output-file',
                    dest='output_file', type="string", default=None,
                    help='Destination filename for reflections')

  # Parse the arguments
  options, args = parser.parse_args()

  # Print help if not enough arguments specified, otherwise call function
  if len(args) < 1:
    parser.print_help()
  else:
    # Initialise the script runner
    runner = ScriptRunner(options, reflections_filename=args[0])

    # Run the script
    runner()
