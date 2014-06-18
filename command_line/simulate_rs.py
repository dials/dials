#!/usr/bin/env python
#
# reciprocal_space.py
#
#  Copyright (C) 2014 Diamond Light Source, James Parkhurst
#
#  This code is distributed under the BSD license, a copy of which is
#  included in the root directory of this package.

from __future__ import division

if __name__ == '__main__':

  from optparse import OptionParser
  from dials.algorithms.simulation.reciprocal_space import Simulator
  from dxtbx.model.experiment.experiment_list import ExperimentListFactory
  from dials.util.command_line import Command
  from math import pi

  usage = "usage: %prog [options] experiment.json"
  parser = OptionParser(usage)

  # Write the datablock to JSON or Pickle
  parser.add_option(
    "-o", "--output",
    dest = "output",
    type = "string", default = "simulated.pickle",
    help = "The output pickle file")

  # Write the datablock to JSON or Pickle
  parser.add_option(
    "-n", "--num",
    dest = "num",
    type = "int", default = 1000,
    help = "The number of reflections (default 1000)")

  # Write the datablock to JSON or Pickle
  parser.add_option(
    "-i", "--intensity",
    dest = "intensity",
    type = "int", default = 1000,
    help = "The intensity (default 1000)")

  # Write the datablock to JSON or Pickle
  parser.add_option(
    "-b", "--background",
    dest = "background",
    type = "string", default = "0",
    help = "The background (default 0)")

  # Write the datablock to JSON or Pickle
  parser.add_option(
    "-r", "--random",
    dest = "random",
    action = "store_true", default=False,
    help = "Random intensities and backgrounds")

  # The profile parameters
  parser.add_option(
    "--sigma_b",
    dest = "sigma_b",
    type = "float", default = None,
    help = "Sigma B")

  parser.add_option(
    "--sigma_m",
    dest = "sigma_m",
    type = "float", default = None,
    help = "Sigma M")

  parser.add_option(
    "--n_sigma",
    dest = "n_sigma",
    type = "float", default = 3,
    help = "Number of standard deviations for profile")

  # Parse the command line arguments
  (options, args) = parser.parse_args()

  # Ensure we have enough arguments
  if len(args) != 1:
    parser.show_help()
    exit(0)

  # Check we have some profile parameters
  if options.sigma_b is None or options.sigma_m is None:
    print 'Need to set profile parameters'
    exit(0)

  print 'Simulating with the following parameters:'
  print ' # Reflections: %d' % options.num
  print ' Intensity: %d' % options.intensity
  print ' Background: %s' % options.background
  print ' Sigma B: %f degrees' % options.sigma_b
  print ' Sigma M: %f degrees' % options.sigma_m
  print ' N Sigma: %f degrees' % options.n_sigma
  print ' Random: %s' % str(options.random)

  # Get the experiments
  experiments = ExperimentListFactory.from_json_file(args[0], check_format=False)
  assert(len(experiments) == 1)
  experiment = experiments[0]

  # Do the simulation
  sigma_b = options.sigma_b * pi / 180
  sigma_m = options.sigma_m * pi / 180
  n_sigma = options.n_sigma
  N = options.num
  I = options.intensity
  B = map(int, options.background.split(","))
  B = B + [0] * (4 - len(B))
  simulate = Simulator(experiment, sigma_b, sigma_m, n_sigma)
  if options.random:
    refl = simulate.with_random_intensity(N, I, B[0], B[1], B[2], B[3])
  else:
    refl = simulate.with_given_intensity(N, I, B[0], B[1], B[2], B[3])

  # Save the reflections to file
  Command.start('Writing reflections to %s' % options.output)
  refl.as_pickle(options.output)
  Command.end('Write reflections to %s' % options.output)
