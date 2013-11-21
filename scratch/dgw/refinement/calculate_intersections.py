#!/usr/bin/env python
#
# calculate_intersections.py
#
#  Copyright (C) (2013) STFC Rutherford Appleton Laboratory, UK.
#
#  Author: David Waterman.
#
#  This code is distributed under the BSD license, a copy of which is
#  included in the root directory of this package.

"""This script calculates and sets the intersections of a set of reflections
with a detector. The reflections must have valid beam_vectors already set!
Reflections whose beam_vectors are (0,0,0) are ignored."""

from __future__ import division

class ScriptRunner(object):
  """Class to run script."""

  def __init__(self, options, sweep_filename, reflections_filename):
    """Setup the script."""

    # Filename data
    self.output_filename = options.output_file
    self.sweep_filename = sweep_filename
    self.reflections_filename = reflections_filename

  def __call__(self):
    """Run the script."""
    from dials.model.serialize import load, dump
    from dials.model.data import ReflectionList
    import cPickle as pickle
    from dials.algorithms.spot_prediction import ray_intersection

    # Load the reflection list
    print 'Loading reflections from {0}'.format(self.reflections_filename)
    rlist = pickle.load(open(self.reflections_filename, 'r'))

    # Try to load the models
    print 'Loading models from {0}'.format(self.sweep_filename)

    sweep = load.sweep(open(self.sweep_filename, 'r'))
    beam = sweep.get_beam()
    wavelength = beam.get_wavelength()
    detector = sweep.get_detector()

    # get the intersections
    observations = ray_intersection(detector, rlist)

    if len(observations) != len(rlist):
      print "WARNING: not all reflections intersect the detector"

      # Why is this? Dump out the unique reflections to explore
      unique = ReflectionList()
      for r in rlist:
        try:
          obs = ray_intersection(detector, r)
        except RuntimeError:
          unique.append(r)

      unique_filename = "unique.pickle"
      print 'Those reflections that do not intersect have been saved' \
            ' to {0}'.format(unique_filename)
      pickle.dump(observations, open(unique_filename, 'wb'),
          pickle.HIGHEST_PROTOCOL)

    # update the centroid positions too
    for r in observations:
      r.centroid_position = r.image_coord_mm + (r.rotation_angle, )

    # Write out reflections
    if self.output_filename is not None:

      print 'Saving reflections to {0}'.format(self.output_filename)
      pickle.dump(observations, open(self.output_filename, 'wb'),
          pickle.HIGHEST_PROTOCOL)

if __name__ == '__main__':

  from optparse import OptionParser

  # Specify the command line options
  usage  = "usage: %prog [options] sweep.json reflections.pickle"

  # Create an option parser
  parser = OptionParser(usage)

  parser.add_option('-o', '--output-file',
                    dest='output_file', type="string", default=None,
                    help='Destination filename for reflections')

  # Parse the arguments
  options, args = parser.parse_args()

  # Print help if not enough arguments specified, otherwise call function
  if len(args) < 2:
    parser.print_help()
  else:
    # Initialise the script runner
    runner = ScriptRunner(options, sweep_filename=args[0],
                          reflections_filename=args[1])

    # Run the script
    runner()
